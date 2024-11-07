import requests
import urllib.request
import shutil
import json
from protocol.taxonomy import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
    
############ Reading files #################

def read_polymers_file(type):
    file_path = f'data/polymers/{type}_records.json'
    with open(file_path, 'r', encoding='utf-8-sig') as file:
        data = json.load(file)
    
    try:
        entities = data[0]["collect(properties(m))"]
        ids = set()
        chains = []
        for obj in entities:
            if obj['parent_rcsb_id'] not in ids:
                ids.add(obj['parent_rcsb_id'])
                chains.append({'parent_id': obj['parent_rcsb_id'], 
                                 'auth_asym_id': obj['auth_asym_id'], 
                                 'tax_id': obj['src_organism_ids'],
                                 'seq': obj['entity_poly_seq_one_letter_code_can']})
             
        return chains   
    except:
        print("Error accessing polymer data.")
        return None
    
def read_polymers_file_filter_by_kingdom(type, kingdom):
    file_path = f'data/polymers/{type}_records.json'
    with open(file_path, 'r', encoding='utf-8-sig') as file:
        data = json.load(file)
    
    try:
        entities = data[0]["collect(properties(m))"]
        ids = set()
        chains = []
        for obj in entities:
            if obj['parent_rcsb_id'] not in ids and TaxId.superkingdom(obj['src_organism_ids'][0]) == kingdom:
                ids.add(obj['parent_rcsb_id'])
                chains.append({'parent_id': obj['parent_rcsb_id'], 
                                 'auth_asym_id': obj['auth_asym_id'], 
                                 'tax_id': obj['src_organism_ids'],
                                 'seq': obj['entity_poly_seq_one_letter_code_can']})
             
        return chains   
    except:
        print("Error accessing polymer data.")
        return None
    
def get_uniprot_seq_from_file(type, rcsb_id):
    if type == "25-28S":
        file_path = f'data/polymers/25S_records.json'
        with open(file_path, 'r', encoding='utf-8-sig') as file:
            data_25 = json.load(file)
        
        file_path = f'data/polymers/28S_records.json'
        with open(file_path, 'r', encoding='utf-8-sig') as file:
            data_28 = json.load(file)
         
        try:
            entities = data_25[0]["collect(properties(m))"]
            for obj in entities:
                if obj['parent_rcsb_id'] == rcsb_id:
                    return obj['entity_poly_seq_one_letter_code_can']
            
            entities = data_28[0]["collect(properties(m))"]
            for obj in entities:
                if obj['parent_rcsb_id'] == rcsb_id:
                    return obj['entity_poly_seq_one_letter_code_can']
            
            return None
        except:
            print("Error accessing polymer data.")
            return None   
        
    else:       
        file_path = f'data/polymers/{type}_records.json'
        with open(file_path, 'r', encoding='utf-8-sig') as file:
            data = json.load(file)
    
    try:
        entities = data[0]["collect(properties(m))"]
        for obj in entities:
            if obj['parent_rcsb_id'] == rcsb_id:
                return obj['entity_poly_seq_one_letter_code_can']
        return None
    except:
        print("Error accessing polymer data.")
        return None
    
def retrieve_taxid(rcsb_id: str):
    file_path = f'data/polymers/uL4_records.json'
    with open(file_path, 'r', encoding='utf-8-sig') as file:
        data = json.load(file)
    
    try:
        entities = data[0]["collect(properties(m))"]
        for obj in entities:
            if obj['parent_rcsb_id'] == rcsb_id:
                return list(obj['src_organism_ids'])[0]
        return None
    except:
        print("Error accessing polymer data.")
        return None
    
def find_kingdom(rcsb_id: str):
    tax_id = retrieve_taxid(rcsb_id)
    
    if tax_id is None:
        print(f"Cannot access rcsb_id {rcsb_id} from dataset.")
        return None
    
    kingdom = TaxId.superkingdom(int(tax_id))
    
    return kingdom

############ RiboXYZ ###################

def get_profile(rcsb_id: str):
    try:
        url = f"https://api.ribosome.xyz/structures/profile?rcsb_id={rcsb_id}"
        response = requests.get(url)
        
        polymers = []
        if response.status_code == 200:
            data = response.json()
            for obj in data['proteins']:
                if obj['assembly_id'] == 0:
                    polymers.append({'polymer': obj['nomenclature'][0], 'auth_asym_id': obj['auth_asym_id'], 'seq': obj['entity_poly_seq_one_letter_code_can']})
                
            return polymers
        else:
            print('Error:', response.status_code)
            return None
    except:
        return None
   
############ PDB downloading #################
 
# Call to get mmcif file
def get_mmcif(parent_id):
    try:
        url = f"https://files.rcsb.org/download/{parent_id}.cif"

        with urllib.request.urlopen(url) as response, open(f'data/mmcif/{parent_id}.cif', 'wb') as out_file:
            shutil.copyfileobj(response, out_file) 
        return True    
    except urllib.error.HTTPError as e:
        print(f"Failed to download {parent_id}.cif: HTTP error {e.code}")
    except urllib.error.URLError as e:
        print(f"Failed to download {parent_id}.cif: URL error {e.reason}")
    except Exception as e:
        print(f"An unexpected error occurred while downloading {parent_id}.cif: {e}")
    return False


def get_sequences_from_PDB(parent_id: str):
    try:
        url = f"https://www.rcsb.org/fasta/entry/{parent_id}"

        with urllib.request.urlopen(url) as response, open(f'data/fasta/parent_files/{parent_id}.fasta', 'wb') as out_file:
            shutil.copyfileobj(response, out_file) 
        return True    
    except urllib.error.HTTPError as e:
        print(f"Failed to download {parent_id} fasta file: HTTP error {e.code}")
    except urllib.error.URLError as e:
        print(f"Failed to download {parent_id} fasta file: URL error {e.reason}")
    except Exception as e:
        print(f"An unexpected error occurred while downloading {parent_id} fasta file: {e}")
    return False

    
############ Working with Fasta ############## 
     
def save_fasta(sequences, type):
    if sequences is None:
        print('No sequences to save.')
        return
    
    with open(f"data/output/fasta/sequences_{type}.fasta", "w") as file:
        for seq in sequences:
            file.write(f">{type}_{seq['parent_id']}_{seq['auth_asym_id']}\n{seq['seq']}\n")
    
            
def save_fasta_by_kingdom(sequences, type, kingdom):
    if sequences is None:
        print('No sequences to save.')
        return
    
    with open(f"data/fasta/sequences_{kingdom}_{type}.fasta", "w") as file:
        for seq in sequences:
            file.write(f">{type}_{seq['parent_id']}_{seq['auth_asym_id']}\n{seq['seq']}\n")
            

def check_fasta_for_rcsb_id(rcsb_id: str, polymer: str, kingdom: str):
    path = f"data/fasta/sequences_{kingdom}_{polymer}.fasta"
    for record in SeqIO.parse(path, "fasta"):
        if rcsb_id in record.id:
            return True
    return False

def add_polymer_to_fasta_list(rcsb_id: str, polymer_type: str, profile: list[dict], kingdom:str):
    polymer = next((item for item in profile if item["polymer"] == polymer_type), None)
    asym_id = polymer['auth_asym_id']
    seq = polymer['seq']
    name = f"{polymer_type}_{rcsb_id}_{asym_id}"
    
    # Create a SeqRecord object
    record = SeqRecord(Seq(seq), id=name, description="")

    # Specify the file path
    file_path = f"data/fasta/sequences_{kingdom}_{polymer_type}.fasta"

    # Parse existing FASTA and add records to list
    records = []
    for seq_rec in SeqIO.parse(file_path, "fasta"):
        records.append(seq_rec)
        
    # Add new record to list  
    records.append(record)
        
    # Rewrite FASTA file including new record
    with open(file_path, "w") as file:
        SeqIO.write(records, file, "fasta")
    
    return f"sequences_{kingdom}_{polymer_type}"