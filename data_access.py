import requests

# API call

def get_proteins(type):
    try:
        url = f"https://api.ribosome.xyz/structures/list_polymers_filtered_by_polymer_class?polymer_class={type}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()

            proteins = [{'parent_id': obj['parent_rcsb_id'], 'seq': obj['entity_poly_seq_one_letter_code']} for obj in data['polymers']]
            print(len(proteins))
            return proteins
        else:
            print('Error:', response.status_code)
            return None
    except:
        return None
    
def save_fasta(sequences, type):

    if sequences is None:
        print('No sequences to save. Exiting.')
        return
    
    with open(f"data/input_sequences_{type}.fasta", "w") as file:
        for seq in sequences:
            file.write(f">{type}_{seq['parent_id']}\n{seq['seq']}\n")
            
