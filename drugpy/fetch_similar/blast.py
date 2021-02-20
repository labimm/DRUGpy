import webbrowser
from ftplib import FTP
from functools import lru_cache

import requests
from pymol import cmd as pm

from ..commons import rscript, PLUGIN_DATA_DIR


@pm.extend
def fetch_similar_blast_update():
    """
    The cluster database needs to be updated before the first use of the
    fetch_similar_blast feature. The database ftp://resources.rcsb.org/sequence/clusters/
    is updated weekly, for new data is prudent to run this command weekly.
    """
    rscb_server = FTP("resources.rcsb.org")
    rscb_server.login()

    rscb_server.cwd("/sequence/clusters/")
    for cluster_fname in [
        "bc-100.out",
        "bc-95.out",
        "bc-90.out",
        "bc-70.out",
        "bc-50.out",
        "bc-40.out",
        "bc-30.out",
    ]:
        with open(PLUGIN_DATA_DIR + "/" + cluster_fname, "wb") as cluster_file:
            rscb_server.retrbinary(
                "RETR " + cluster_fname, cluster_file.write, blocksize=262144
            )


def find_similar_chain_ids(chain_id, threshold):
    """Fetch structure similar chain ids from RCSB PDB.
    chain_id - reference chain id
    threshold - similarity threshold
    """

    cluster_fname = f"bc-{threshold}.out"
    with open(PLUGIN_DATA_DIR + "/" + cluster_fname) as cluster_file:
        for cluster in cluster_file:
            if chain_id.upper() in cluster.upper():
                break
        else:
            return []

        sim_chain_ids = []
        for chain_id in cluster.split():
            pdb, chain = chain_id.split("_")
            sim_chain_ids.append((pdb.upper(), chain.upper()))
        return sim_chain_ids


@lru_cache()
def get_resolution(pdb_id):
    """
    Get the resolution for a PDB id, or None case it doesn't have.
    """
    ret = requests.post(
        "https://data.rcsb.org/graphql",
        json={
            "query": f"""
    {{
        entry(entry_id: "{pdb_id}") {{
            pdbx_vrpt_summary {{
                PDB_resolution
            }}
      }}
    }}
    """
        },
    )
    data = ret.json()
    resol = data["data"]["entry"]["pdbx_vrpt_summary"]["PDB_resolution"]
    return resol


def plot_hierarquical_cluster(chain_ids):
    chain_ids = ",".join(map(repr, chain_ids))
    out, success = rscript(
        f"""
        library(bio3d)
        ids <- c({chain_ids})
        files <- get.pdb(ids, path="pdbs", split=TRUE)
        pdbs <- pdbaln(files)
        cores <- core.find(pdbs)
        xyz <- pdbfit(pdbs, inds=cores)
        pc <- pca(xyz, rm.gaps=TRUE)
        d <- dist(pc$z[, 1:2])
        hc <- hclust(d)
        hclustplot(hc, k=1, labels = ids)
    """
    )
    webbrowser.open("Rplots.pdf")
    print(out)


@pm.extend
def fetch_similar_blast(
    chain_id,
    similarity=95,
    ligand=None,
    dist=5,
    compounds="organic or inorganic",
    prosthetic_groups="HEM FAD NAP NDP ADP FMN",
    max_resolution=None,
    max_structures=50,
):
    """
    Fetch sequence similar structures from RCSB PDB and optionally keep only
    apo structures. Apo are evaluated respective to a choosen ligand on the
    reference chain.

    On the first use update the database with the command `update_cluster_data`.
    Update the database weekly.

    OPTIONS
        chain_id        Reference structure chain id.
        similarity      Sequence similarity threshold (one of the available
                        from RCSB PDB).
        ligand          Reference ligand PDB id.
        dist            Distance cut-off around reference ligand for apo
                        evaluation.
        compounds       Selection that should be considered ligands upon apo
                        computation. Only used when ligand is given.
        prothestic_groups   List of ligands to be ignored when evaluating apo.
        max_resolution  Fetch only X-ray structures with up to such
                        resolution.
        max_structures  Fetch at most n structures. 0 for all structures.
    EXAMPLES
        fetch_similar_blast 2XY9_A, 100
        fetch_similar_blast 2XY9_A, 95, 3ES, 3, organic
        fetch_similar_blast 6Y2F_A, max_structures=0
    SEE ALSO
        update_cluster_data
        fetch_similar_shape3d
    """
    chain_id = chain_id.upper()
    max_structures = int(max_structures)
    pm.fetch(chain_id, chain_id)

    sims = []
    similars = find_similar_chain_ids(chain_id, similarity)
    cont = 0
    for sim_pdb, sim_chain in similars:

        if max_structures != 0 and cont >= max_structures:
            break

        sim_obj = f"{sim_pdb}_{sim_chain}"
        if sim_obj.upper() == chain_id.upper():
            continue

        pm.fetch(sim_obj, **{"async": 0})
        pm.align(sim_obj, chain_id)

        # Check the resolution
        resol = None
        if max_resolution:
            resol = get_resolution(sim_pdb)
            if not resol or resol > max_resolution:
                pm.delete(sim_obj)
                continue

        # Check nearby ligands
        if ligand:
            model = pm.get_model(
                f"({sim_obj} and ({compounds}))"
                f" within {dist} of"
                f"({chain_id} and (resn {ligand}))"
            )
            resns = set(a.resn for a in model.atom)

            is_apo = True
            for resn in resns:
                if resn not in prosthetic_groups.split():
                    is_apo = False
                    break

            if not is_apo:
                pm.delete(sim_obj)
                continue

        cont += 1
        sims.append((sim_obj, sim_chain, sim_pdb, resol))

    plot_hierarquical_cluster([chain_id] + [s[0] for s in sims])
    return sims
