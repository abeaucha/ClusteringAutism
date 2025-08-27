#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# fetch_reactome_hierarchy.py
# Author: Antoine Beauchamp
# Created: August 27th, 2025

# Packages -------------------------------------------------------------------

import os
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import time
import pandas as pd


# Environment variables ------------------------------------------------------

PROJECTPATH = os.getenv("PROJECTPATH")


# Main -----------------------------------------------------------------------

if __name__ == '__main__':

    session = requests.Session()
    retries = Retry(total = 5, backoff_factor = 1,
                    status_forcelist = [500, 502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries = retries))

    # Get pathways at the top level of the Reactome database
    pathways_root_url = 'https://reactome.org/ContentService/data/pathways/top/9606'
    pathways_root_resp = requests.get(pathways_root_url).json()

    # Combine root names and IDs into a data frame
    df_pathways_root = pd.DataFrame(
        dict(Name = [x['displayName'] for x in pathways_root_resp],
             ID = [x['stId'] for x in pathways_root_resp])
    )

    # Iterate over pathway trees
    list_pathways_all = []
    for i, row in df_pathways_root.iterrows():

        print("Fetching hierarchy for pathway root: {}".format(row['Name']))

        # For each root pathway, get all the events that fall within that tree
        pathway_hierarchy_url = 'https://reactome.org/ContentService/data/pathway/{}/containedEvents/'.format(row['ID'])
        try:
            pathway_hierarchy_resp = session.get(pathway_hierarchy_url, timeout = 10)
            pathway_hierarchy_json = pathway_hierarchy_resp.json()
        except Exception as e:
            print(f"Failed to fetch hierarchy for pathway {row['Name']}: {e}")
            time.sleep(5)

        df_pathway_hierarchy = pd.DataFrame(
            dict(Name = [x['displayName'] for x in pathway_hierarchy_json],
                 ID = [x['stId'] for x in pathway_hierarchy_json],
                 SchemaClass = [x['schemaClass'] for x in pathway_hierarchy_json])
        )


        print("\tFetching parents of nodes in hierarchy...")

        # For each root pathway, get the parent pathway of each event
        pathway_parents_url = '{}eventOf'.format(pathway_hierarchy_url)
        try:
            pathway_parents_resp = session.get(pathway_parents_url, timeout = 10)
        except Exception as e:
            print(f"Failed to fetch parents of nodes in hierarchy: {e}")
            time.sleep(5)

        pathway_parents = (pathway_parents_resp.text
                           .replace("[", "").replace("[", "")
                           .split("\n, "))
        pathway_parents = [pathway.split("\n") for pathway in pathway_parents]
        pathway_parents = [pathway[0] for pathway in pathway_parents]
        pathway_parents = [pathway.split("\t") for pathway in pathway_parents]
        df_pathway_parents = pd.DataFrame(
            dict(ParentID = [x[0] for x in pathway_parents],
                 Parent = [x[1] for x in pathway_parents])
        )

        print("\tJoining nodes with parent nodes...")

        # Concatenate node information with parent node information
        df_pathways = pd.concat([df_pathway_hierarchy, df_pathway_parents], axis = 1)
        df_pathways['Root'] = row['Name']
        df_pathways = df_pathways.loc[df_pathways['SchemaClass'] == 'Pathway'].copy()
        df_pathways = df_pathways.drop("SchemaClass", axis = 1)

        # Create a row for the hierarchy root
        df_pathway_root = pd.DataFrame([
            dict(Name = row['Name'], ID = row['ID'],
                 ParentID = None, Parent = None,
                 Root = row['Name'])
        ])

        # Concatenate the root with the rest of the hierarchy
        df_pathways = pd.concat([df_pathway_root, df_pathways], axis = 0)
        df_pathways = df_pathways.reset_index(drop = True)

        # Append hierarchy to list
        list_pathways_all.append(df_pathways)

# Concatenate all pathways into one data frame
df_pathways_all = pd.concat(list_pathways_all, axis = 0).reset_index(drop = True)

# Export pathway hierarchy
outfile = "reactome_hierarchy.csv"
outdir = os.path.join(PROJECTPATH, "data/enrichment/")
outfile = os.path.join(outdir, outfile)
df_pathways_all.to_csv(outfile, index = False)