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

# Create an RDataFrame with all files in the dataset
# Replace "events" with the actual tree name
rdf = ROOT.RDataFrame("events", file_paths)

# Create histogram of EventHeader.eventNumber
h_eventNumber = rdf.Histo1D(
    ("h_eventNumber", "Event Number Distribution;Event Number;Count", 100, 0, 100),
    "EventHeader.eventNumber"
)

# Draw and save the histogram
canvas = ROOT.TCanvas("c1", "Event Number Distribution", 800, 600)
h_eventNumber.Draw()
canvas.SaveAs("event_number_distribution.png")

# Print total entries
print(f"Total entries in dataset: {rdf.Count().GetValue()}")
print("Saved histogram to event_number_distribution.png")
