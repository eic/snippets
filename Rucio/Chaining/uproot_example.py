from rucio.client import Client
import uproot
import matplotlib.pyplot as plt
import numpy as np

# Initialize Rucio client
client = Client()

# Define the dataset DID
dataset_did = "epic:/RECO/26.02.0/epic_craterlake/SINGLE/e+/500MeV/3to50deg"
scope, name = dataset_did.split(':', 1)

# Get the list of files in the dataset
files = list(client.list_files(scope, name))
dids = [{'scope': f['scope'], 'name': f['name']} for f in files]

# Get one replica PFN for each file in the dataset
file_paths = [
    next(iter(replica['pfns']))  # Get first PFN URL (dict keys are the URLs)
    for replica in client.list_replicas(dids, rse_expression='isopenaccess=true')
    if replica['pfns']
]

# Collect EventHeader.eventNumber from all files
event_numbers = []
for file_path in file_paths:
    with uproot.open(file_path) as f:
        tree = f["events"]  # Replace "events" with the actual tree name
        event_numbers.extend(tree["EventHeader.eventNumber"].array())

# Create histogram
plt.hist(event_numbers, bins=50)
plt.xlabel("Event Number")
plt.ylabel("Count")
plt.title("EventHeader.eventNumber Distribution")
plt.savefig("event_number_distribution.png")
print("Saved histogram to event_number_distribution.png")
