from rucio.client import Client
import ROOT

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

# Create a TChain to process all files as a single dataset
# Replace "events" with the actual tree name in your files
chain = ROOT.TChain("events")

for file_path in file_paths:
    chain.Add(file_path)

# Create histogram of EventHeader.eventNumber
canvas = ROOT.TCanvas("c1", "Event Number Distribution", 800, 600)
chain.Draw("EventHeader.eventNumber>>h_eventNumber(100)")
canvas.SaveAs("event_number_distribution.png")

# Print total entries
print(f"Total entries in dataset: {chain.GetEntries()}")
print("Saved histogram to event_number_distribution.png")
