import subprocess

# List of 50 rcsb_id's
rcsb_ids = ['1S72', '8HL4', '4ADX', '3I55', '3CC4', '1VQ8', 
            '1M1K', '3G71', '1VQN', '1VQ5', '1N8R', '2QEX', 
            '1YJN', '8HL2', '1K8A', '1Q81', '3CCJ', '8HL5', 
            '1KC8', '3G6E', '3CME', '2OTL', '8HL1', '1Q86', 
            '3CC2', '6SKF', '1VQ7', '3CMA', '1Q7Y', '6TH6', 
            '1VQK', '3CCE', '3I56', '1KD1', '3CD6', '1K73', 
            '8HL3', '3CCS', '4V4N', '1M90', '3CXC', '1VQ4', 
            '8HKV', '8HKU', '1KQS', '2OTJ', '3CCQ', '1FFK', 
            '3CCR', '3OW2', '1JJ2']

# Iterate over each rcsb_id and run main.py
for rcsb_id in rcsb_ids:
    # Run main.py with the current rcsb_id as a parameter
    result = subprocess.run(
        ["python3", "-m", "protocol.main", rcsb_id],
        capture_output=True,
        text=True
    )
    
    # Print the output and error (if any)
    print(f"Output for RCSB ID {rcsb_id}:\n{result.stdout}")
    if result.stderr:
        print(f"Error for RCSB ID {rcsb_id}:\n{result.stderr}")
