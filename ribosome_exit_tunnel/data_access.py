import requests
import urllib.request
import shutil
import json
import time

# API call to get list of all proteins of type 'type'
def get_proteins(type):
    try:
        url = f"https://api.ribosome.xyz/structures/list_polymers_filtered_by_polymer_class?polymer_class={type}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            
            proteins = []
            ids = set()
            for obj in data['polymers']:
                if obj['parent_rcsb_id'] not in ids:
                    ids.add(obj['parent_rcsb_id'])
                    proteins.append({'parent_id': obj['parent_rcsb_id'], 'auth_asym_id': obj['auth_asym_id'], 'seq': obj['entity_poly_seq_one_letter_code']})
                
            return proteins
        else:
            print('Error:', response.status_code)
            return None
    except:
        return None
    
def read_polymers_file(type):
    file_path = f'data/polymers/{type}_records.json'
    with open(file_path, 'r', encoding='utf-8-sig') as file:
        data = json.load(file)
    
    try:
        entities = data[0]["collect(properties(m))"]
        ids = set()
        proteins = []
        for obj in entities:
            if obj['parent_rcsb_id'] not in ids:
                ids.add(obj['parent_rcsb_id'])
                proteins.append({'parent_id': obj['parent_rcsb_id'], 'auth_asym_id': obj['auth_asym_id'], 'seq': obj['entity_poly_seq_one_letter_code']})
             
        return proteins   
    except:
        print("Error accessing polymer data.")
        return None
    
    return
    
# API call to get mmcif file
def get_mmcif(parent_id):
    try:
        url = f"https://files.rcsb.org/download/{parent_id}.cif"

        with urllib.request.urlopen(url) as response, open(f'data/mmcif/{parent_id}.cif', 'wb') as out_file:
            shutil.copyfileobj(response, out_file)     
    except:
        return None
    
   
def save_fasta(sequences, type):
    if sequences is None:
        print('No sequences to save.')
        return
    
    with open(f"data/output/input_sequences_{type}.fasta", "w") as file:
        for seq in sequences:
            file.write(f">{type}_{seq['parent_id']}_{seq['auth_asym_id']}\n{seq['seq']}\n")
            