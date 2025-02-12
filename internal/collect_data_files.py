import csv
from pathlib import Path

def collect_landmarks_into_file(kingdom: str = None):

    landmarks = []
    
    if kingdom is None:
        path1 = "data/output/conserved/alignment_close_conserved.csv"
        kingdom = "universal"
    else:
        path1 = f"data/output/conserved/alignment_close_conserved_{kingdom}.csv"

    with open(path1, mode='r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            landmarks.append(f"{row['chain']}-x")
            landmarks.append(f"{row['chain']}-y")
            landmarks.append(f"{row['chain']}-z")

    with open(f"data/output/landmarks/{kingdom}.csv", mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames= ['rcsb_id'] + landmarks)
        
        if file.tell() == 0:
            writer.writeheader()
        
        folder_path = Path(f"data/output/landmarks/{kingdom}")
        
        for file in folder_path.iterdir():
            with open(file, mode='r', newline='') as file_i:
                reader = csv.DictReader(file_i)
                row_to_write = {}
                for row in reader:
                    row_to_write['rcsb_id'] = row['parent_id']
                    row_to_write[f"{row['landmark']}-x"] = int(row['x'])
                    row_to_write[f"{row['landmark']}-y"] = int(row['y'])
                    row_to_write[f"{row['landmark']}-z"] = int(row['z'])
                writer.writerow(row_to_write)