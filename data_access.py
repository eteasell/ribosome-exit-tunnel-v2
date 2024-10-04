import requests
import urllib.request
import shutil

# API call to get list of all proteins of type 'type'
def get_proteins(type):
    try:
        url = f"https://api.ribosome.xyz/structures/list_polymers_filtered_by_polymer_class?polymer_class={type}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()

            proteins = [{'parent_id': obj['parent_rcsb_id'], 'auth_asym_id': obj['auth_asym_id'], 'seq': obj['entity_poly_seq_one_letter_code']} 
                        for obj in data['polymers']]
            return proteins
        else:
            print('Error:', response.status_code)
            return None
    except:
        return None
    
# API call to get mmcif file for given protein chain
def get_mmcif(parent_id, chain):
    try:
        url = f"https://api.ribosome.xyz/mmcif/polymer?rcsb_id={parent_id}&auth_asym_id={chain}"

        with urllib.request.urlopen(url) as response, open(f'data/mmcif/{parent_id}_{chain}.cif', 'wb') as out_file:
            shutil.copyfileobj(response, out_file)     
    except:
        return None
    
def save_fasta(sequences, type):
    if sequences is None:
        print('No sequences to save.')
        return
    
    with open(f"data/input_sequences_{type}.fasta", "w") as file:
        for seq in sequences:
            file.write(f">{type}_{seq['parent_id']}_{seq['auth_asym_id']}\n{seq['seq']}\n")
            