import pandas as pd
import uuid
import sys
import shlex  # Import shlex for quoting

if len(sys.argv) != 3:
    print("Usage: python generate_batch_commands.py <path_to_csv> <path_to_new_csv>")
    sys.exit(1)

csv_path = sys.argv[1]
csv_new_path = sys.argv[2]

# Load the sessions CSV
sessions_df = pd.read_csv(csv_path)

# Generate a UUID for each row
sessions_df['UUID'] = [str(uuid.uuid4()) for _ in range(len(sessions_df))]

# Add UUID as the first column
cols = sessions_df.columns.tolist()
cols = cols[-1:] + cols[:-1]
sessions_df = sessions_df[cols]

# Generate sbatch commands and include the UUID in the command
for index, row in sessions_df.iterrows():
    session_details = row.to_json()
    session_details_quoted = shlex.quote(session_details)  # Properly quote the JSON string
    sbatch_command = f"sbatch --job-name=job_{row['UUID']} ece_run.sh {session_details_quoted}"
    print(sbatch_command)

# Save the updated DataFrame to a new CSV file or overwrite the existing one
sessions_df.to_csv(csv_new_path, index=False)
