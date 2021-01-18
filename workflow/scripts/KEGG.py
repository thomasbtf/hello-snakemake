import pandas as pd
import urllib.parse
import urllib.request
import time
from typing import List

def chunks(lst:list, n:int):
    """Yield successive n-sized chunks from list.

    Args:
        lst (list): Orignal list, that shall be chunked
        n (int): Size of chunks

    Yields:
        list: n-sized chunks from lst
    """

    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def map_db_identifiers(from_id:str, to_id:str, query:List[str], verbose:bool = True) -> pd.DataFrame:
    """Querys the UniProt Retrieve/ID mapping servies to map database identifiers. Converts identifiers 
    which are of a different type to UniProt identifiers or vice versa and returns the identifier lists. 
    See https://www.uniprot.org/help/api_idmapping for further information.
    Alternativly load the mapping data directly from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping (~19.5 GB).
    Very large mapping requests (>50,000 identifiers) are likely to fail. Please do verify that your list does not contain any duplicates.
    Very large mapping requests are split it into smaller chunks (<20,000 identifiers). 

    Args:
        from_id (str): Original identifier
        to_id (str): Target identifier
        query (List[str]): Identifiers to be mapped
        verbose (bool, optional): Verbose mode. Provides additional details. Defaults to True.

    Returns:
        pd.DataFrame: [description]
    """

    url = 'https://www.uniprot.org/uploadlists/'
    return_df = pd.DataFrame()
    chunk_size = 20000

    # chunk list if to large
    iterationCount = 1
    len_gen = len(list(chunks(query, chunk_size)))
    for chunk in chunks(query, chunk_size):
        chunk_downloaded = False

        params = {
            'from': from_id,
            'to': to_id,
            'format': 'tab',
            'query': " ".join(chunk)
        }
       
        # encode parameters
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')

        # # create request object and read it
        if verbose:
            print('Querying chunk %s/%s' %(iterationCount, len_gen))
        request = urllib.request.Request(url, data)

        # Try to download chuck. If not possible, repeat.
        while not chunk_downloaded:
            try:
                with urllib.request.urlopen(request) as req:
                    response = req.read().decode('utf-8')
                chunk_downloaded = True
            except Exception as e:
                if verbose:
                    print('Error querying chunk %s/%s - trying again'%(iterationCount, len_gen))
                    print(e)
                    time.sleep(10)
                pass

        # store response in dataframe
        temp_df = pd.DataFrame([x.split('\t') for x in response.split('\n')])
        temp_df.rename(columns=temp_df.iloc[0], inplace=True)
        temp_df.drop(temp_df.index[0], inplace=True)
        temp_df.drop(temp_df.index[-1], inplace=True)
        return_df = return_df.append(temp_df, ignore_index=True)
        iterationCount += 1

    return return_df.drop_duplicates()

def get_kegg_annotations(data_path, out_path):
    diamond_headers = [
    'Query accession', # the accession of the sequence that was the search query against the database, as specified in the input FASTA file after the > character until the first blank.
    'Target accession', # the accession of the target database sequence (also called subject) that the query was aligned against.
    'Sequence identity', # The percentage of identical amino acid residues that were aligned against each other in the local alignment.
    'Length', # The total length of the local alignment, which including matching and mismatching positions of query and subject, as well as gap positions in the query and subject.
    'Mismatches', # The number of non-identical amino acid residues aligned against each other.
    'Gap openings', # The number of gap openings.
    'Query start', # The starting coordinate of the local alignment in the query (1-based).
    'Query end', # he ending coordinate of the local alignment in the query (1-based).
    'Target start', # The starting coordinate of the local alignment in the target (1-based).
    'Target end', # The ending coordinate of the local alignment in the target (1-based).
    'E-value', # The expected value of the hit quantifies the number of alignments of similar or better quality that you expect to find searching this query against a database of random sequences the same size as the actual target database. This number is most useful for measuring the significance of a hit. By default, DIAMOND will report all alignments with e-value < 0.001, meaning that a hit of this quality will be found by chance on average once per 1,000 queries.
    'Bit score', # The bit score is a scoring matrix independent measure of the (local) similarity of the two aligned sequences, with higher numbers meaning more similar. It is always >= 0 for local Smith Waterman alignments.
    ]   

    # Load identified proteins
    data_df = pd.read_csv(data_path, names = diamond_headers)
    data_df['Target'] = data_df['Target accession'].str.extract(r'(.*(?=\.))')

    # Get KEGG IDs with UniProt Retrieve/ID mapping
    query_ids = data_df['Target'].unique().tolist()#[:100]
    kegg_ids = map_db_identifiers('ACC+ID', 'KEGG_ID', query_ids)

    # Agg. KEEG IDs into lists
    grouped_eco = kegg_ids[kegg_ids['To'].str.contains('eco:')]
    grouped_eco = grouped_eco.groupby('From')['To'].apply(list)
    grouped_keg = kegg_ids.groupby('From')['To'].apply(list)

    grouped_eco = pd.DataFrame(grouped_eco.reset_index()).rename(columns={'To':'KEGG eco tag'})
    grouped_keg = pd.DataFrame(grouped_keg.reset_index()).rename(columns={'To':'KEGG all tags'})

    # Merge with identified proteins
    data_df = data_df.merge(grouped_eco, how='left', left_on='Target', right_on='From', validate='m:1')
    data_df = data_df.merge(grouped_keg, how='left', left_on='Target', right_on='From', validate='m:1')
    data_df.drop(columns=['From_x', 'From_y'], inplace=True)

    data_df.to_csv(out_path, index=False)

get_kegg_annotations(snakemake.input[0], snakemake.output[0])