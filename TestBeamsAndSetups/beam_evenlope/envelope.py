import math
import ROOT
import time #optional to make plots not appear behind terminal

# --- Generation and Assignment ---
class EventSimulator:

    # Variable Dictionary
    # Additional particles can be added here, and changing the vertex and momentum constraints can be done here
    def __init__(self, output_file, n_events=10000):
        self.output_file = output_file
        self.n_events = n_events
        self.min_mom = 3.0 # Lower limit of momentum - 3 G/eV
        self.max_mom = 10.0 # Upper limit of momentum - 10 G/eV
        self.square_side = 300.0 # Limits vertex region
        self.rng = ROOT.TRandom3(0)
        self.PARTICLE_DATA = {
            "1": {"name": "pion", "pdg_id": 211, "mass": 0.13957039},
            "2": {"name": "kaon", "pdg_id": 321, "mass": 0.493677},
            "3": {"name": "proton", "pdg_id": 2212, "mass": 0.93827208816},
            "5": {"name": "Deuteron", "pdg_id": 1000010020, "mass": 1.87561294257},
            "6": {"name": "Alpha Particle", "pdg_id": 1000020040, "mass": 3.7273794066},


        }

    # Prompt User for Particle Species
    def get_particle_choice(self):
        while True:
            print("\nOptions:")
            print("1. Pion")
            print("2. Kaon")
            print("3. Proton")
            print("4. Mixed Generation")
            print("5. Deuteron (D2 nucleus)")
            print("6. Alpha (He4 nucleus)")
            choice = input("Enter your choice: ")
            if choice in self.PARTICLE_DATA:
                return self.PARTICLE_DATA[choice]
            if choice == "4":
                return "mixed"
            print("Invalid input. Please enter 1, 2, or 3.")

    # Randomly Generates x,y,z momentum in a given range
    def generate_momentum(self):
        while True:
            p = self.rng.Gaus(6.5, 1.0)
            if self.min_mom <= p <= self.max_mom:
                break
        theta = self.rng.Uniform(0, 0.030) #30milliradians
        phi = self.rng.Uniform(0, 2 * math.pi)
        px = p * math.sin(theta) * math.cos(phi)
        py = p * math.sin(theta) * math.sin(phi)
        pz = p * math.cos(theta)
        return px, py, pz, p

    # Randomly Generates x,y vertex coordinates in a given range, this range can be editted via self.square_side =... in the variable dictionary
    def generate_vertex(self):
        while True:
            vx = self.rng.Gaus(0, 50)
            vy = self.rng.Gaus(0, 50)
            if (abs(vx) <= self.square_side / 2) and (abs(vy) <= self.square_side / 2):
                break
        vz = 0.0
        return vx, vy, vz

    # Event details for hepmc output file
    def write_event(self, f, event_id, px, py, pz, E, vx, vy, vz, pdg_id):
         event_str = (
             f"E {event_id} 0 0 0 0 0\n"
             "F GenEvent 3 0\n"
             "F Units 0 0\n"
             "F MomentumUnit GEV\n"
             "F LengthUnit MM\n"
             f"F GenParticle 1 1 {pdg_id} 1 0 0 0 0 {px:.6f} {py:.6f} {pz:.6f} {E:.6f}\n"
             f"F GenVertex 1 -1 -1 -1 {vx:.3f} {vy:.3f} {vz:.3f} 0.0\n"
             "F GenParticle 1\n"
             "F GenVertex 1\n"
          )
         f.write(event_str)

    # New write algorithm tested
    # Particle definition: P <barcode> <status> <PDG ID> <px> <py> <pz> <E> <m> <prod_vertex_barcode>
    # Vertex definition: V <barcode> <status> [<incoming particles>] <x> <y> <z> <t>

    def write_event_new(self, f, event_id, px, py, pz, E, mass, vx, vy, vz, pdg_id):
         event_str = (
             f"E {event_id} 1 2\n"
             "U GEV MM\n"
             f"P 1 0 {pdg_id} {px:.6f} {py:.6f} {pz:.6f} {E:.6f} {mass:.6f} 1\n"
             f"V -1 0 [1] {vx:.3f} {vy:.3f} {vz:.3f} 0.0\n"
             f"P 2 -1 {pdg_id} {px:.6f} {py:.6f} {pz:.6f} {E:.6f} {mass:.6f} 1\n"
          )
         f.write(event_str)

    def generate_events(self):
         choice_result = self.get_particle_choice()

         with open(self.output_file, "w") as f:
             # Output the HepMC header
             f.write("HepMC::Version 3.02.02\n")
             f.write("HepMC::Asciiv3-START_EVENT_LISTING\n")

             if choice_result != "mixed":
                 particle = choice_result
                 mass = particle["mass"]
                 pdg_id = particle["pdg_id"]
                 for event_id in range(self.n_events):
                     px, py, pz, p_mag = self.generate_momentum()
                     vx, vy, vz = self.generate_vertex()
                     E = math.sqrt(p_mag**2 + mass**2)
                     self.write_event_new(f, event_id, px, py, pz, E, mass, vx, vy, vz, pdg_id)

            #Mixed generation
             else:
                # Define particles for the mix
                 proton = self.PARTICLE_DATA["3"]
                 pion = self.PARTICLE_DATA["1"]
                 kaon = self.PARTICLE_DATA["2"]
 
                # Define cumulative probabilities for 5:3:1 ratio, edit this ratio as needed.
                # Total parts = 5 + 3 + 1 = 9
                 prob_proton = 5.0 / 9.0
                 prob_pion_cumulative = prob_proton + (3.0 / 9.0) # Proton + Pion
 
                 for event_id in range(self.n_events):
                     # For each event, randomly select a particle based on the ratio
                     rand_val = self.rng.Uniform(0, 1)
 
                     if rand_val < prob_proton:
                         particle = proton
                     elif rand_val < prob_pion_cumulative:
                         particle = pion
                     else:
                         particle = kaon
 
                     mass = particle["mass"]
                     pdg_id = particle["pdg_id"]
 
                     px, py, pz, p_mag = self.generate_momentum()
                     vx, vy, vz = self.generate_vertex()
                     E = math.sqrt(p_mag**2 + mass**2)
                     self.write_event_new(f, event_id, px, py, pz, E, mass, vx, vy, vz, pdg_id)
 
            # Output the HepMC tail
             f.write("HepMC::Asciiv3-END_EVENT_LISTING\n")
         print(f"Wrote {self.n_events} events to {self.output_file}")
 



# Plotting and Fitting
class PlotManager:
    def __init__(self, input_file):
        self.input_file = input_file
        self.hist_vertex_raw = ROOT.TH2F("vertex", "Vertex Distribution;X (mm);Y (mm)", 100, -200, 200, 100, -200, 200)
        self.hist_momentum_raw = ROOT.TH1F("momentum", "Momentum Magnitude;|p| (GeV);Events", 100, 0, 12)
        self.c_vertex = None
        self.c_momentum = None
        self.plotted_hist_vertex = None
        self.plotted_hist_momentum = None
        self.fit_vertex = None
        self.fit_momentum = None
        self.load_data()

# Loads vertex and momentum data by parsing hepmc file
    def load_data(self):
        with open(self.input_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                    
                # Format: V barcode status [incoming] x y z t
                if parts[0] == 'V' and len(parts) >= 8:
                    try:
                        vx = float(parts[4])
                        vy = float(parts[5])
                        self.hist_vertex_raw.Fill(vx, vy)
                    except (ValueError, IndexError):
                        pass

                # Format: P barcode status pdg_id px py pz E m prod_vtx
                elif parts[0] == 'P' and len(parts) >= 8:
                    try:
                        px = float(parts[4])
                        py = float(parts[5])
                        pz = float(parts[6])
                        p_mag = math.sqrt(px**2 + py**2 + pz**2)
                        self.hist_momentum_raw.Fill(p_mag)
                    except (ValueError, IndexError):
                        pass


    def create_canvases(self, vertex_canvas_name, momentum_canvas_name):
        if self.c_vertex:
            self.c_vertex.Close()
            self.c_vertex = None
        if self.c_momentum:
            self.c_momentum.Close()
            self.c_momentum = None

        canvas_width = 800
        canvas_height = 600
        padding = 30
        momentum_x_pos = 50
        momentum_y_pos = 50
        vertex_x_pos = momentum_x_pos
        vertex_y_pos = momentum_y_pos + canvas_height + padding

        self.c_momentum = ROOT.TCanvas(momentum_canvas_name, "Momentum Plot",
                                       momentum_x_pos, momentum_y_pos,
                                       canvas_width, canvas_height)

        self.c_vertex = ROOT.TCanvas(vertex_canvas_name, "Vertex Plot",
                                     vertex_x_pos, vertex_y_pos,
                                     canvas_width, canvas_height)


    # Plots vertex and momentum histograms
    def plot_histograms(self, fit=False):

         if self.plotted_hist_vertex:
             self.plotted_hist_vertex.Delete()
             self.plotted_hist_vertex = None
         if self.plotted_hist_momentum:
             self.plotted_hist_momentum.Delete()
             self.plotted_hist_momentum = None
         if self.fit_vertex:
             self.fit_vertex.Delete()
             self.fit_vertex = None
         if self.fit_momentum:
             self.fit_momentum.Delete()
             self.fit_momentum = None


         vertex_canvas_name = f"c_vertex_{ROOT.gRandom.Integer(1000000)}"
         momentum_canvas_name = f"c_momentum_{ROOT.gRandom.Integer(1000000)}"

         self.create_canvases(vertex_canvas_name, momentum_canvas_name)


         self.plotted_hist_vertex = self.hist_vertex_raw.Clone(f"vertex_clone_{ROOT.gRandom.Integer(1000000)}")
         self.plotted_hist_momentum = self.hist_momentum_raw.Clone(f"momentum_clone_{ROOT.gRandom.Integer(1000000)}")




            #plot styling
         from ROOT import gStyle
         gStyle.SetOptFit(0)
         gStyle.SetOptStat(0)
         self.plotted_hist_momentum.GetXaxis().CenterTitle()
         self.plotted_hist_momentum.GetYaxis().CenterTitle()
         self.plotted_hist_vertex.GetXaxis().CenterTitle()
         self.plotted_hist_vertex.GetYaxis().CenterTitle()
         self.plotted_hist_vertex.GetYaxis().SetTitleOffset(0.9)




         # 6. Vertex canvas + Fitting
         self.c_vertex.cd()
         self.c_vertex.Clear() # Clear previous draws on the canvas
         self.plotted_hist_vertex.Draw("COLZ")
         self.c_vertex.Update()
         if fit:
             v_fit_name = f"fit2d_{ROOT.gRandom.Integer(100000)}" # Unique name for TF2
             self.fit_vertex = ROOT.TF2(v_fit_name, "[0]*exp(-0.5*((x/[1])**2 + (y/[2])**2))", -150, 150, -150, 150)
             self.fit_vertex.SetParameters(1, 50, 50)
             self.plotted_hist_vertex.Fit(self.fit_vertex, "RS")
             chi2 = self.fit_vertex.GetChisquare()
             ndf = self.fit_vertex.GetNDF()
             reduced_chi2 = chi2 / ndf if ndf > 0 else float('inf')
             print(f"Vertex Fit Reduced χ² = {reduced_chi2:.3f}")
             self.c_vertex.Update()
             self.c_vertex.RaiseWindow() # Optional: Uncomment if you want windows to pop to front
             time.sleep(0.1) # Optional: Uncomment for a brief pause after plot update

         # 7. Momentum canvas + Fitting
         self.c_momentum.cd()
         self.c_momentum.Clear() # Clear previous draws on the canvas
         self.plotted_hist_momentum.Draw()
         self.c_momentum.Update()
         if fit:
             m_fit_name = f"fit1d_{ROOT.gRandom.Integer(100000)}" # Unique name for TF1
             self.fit_momentum = ROOT.TF1(m_fit_name, "gaus", 0, 12)
             self.plotted_hist_momentum.Fit(self.fit_momentum, "RS")
             chi2 = self.fit_momentum.GetChisquare()
             ndf = self.fit_momentum.GetNDF()
             print(f"Momentum Fit Reduced χ²/NDF = {chi2 / ndf:.3f}")
             self.c_momentum.Update()
             self.c_momentum.RaiseWindow() # Optional: Uncomment if you want windows to pop to front
             time.sleep(0.1) # Optional: Uncomment for a brief pause after plot update





# Allows user to view unfitted or fitted plots continuously
    def interactive_plot(self):
        while True:
            print("\nPlot Options:")
            print("1. Unfitted Distributions")
            print("2. Fitted Distributions")
            print("0. Exit")
            choice = input("Enter your choice: ")

            #this loop is to make sure if a plot is X'd out, and a replot is attempted, a plot will be generated
            if choice == "1":
                self.plot_histograms(fit=False)
            elif choice == "2":
                self.plot_histograms(fit=True)
            elif choice == "0":
                print("Exiting...")
                if self.plotted_hist_vertex:
                    self.plotted_hist_vertex.Delete()
                    self.plotted_hist_vertex = None
                if self.plotted_hist_momentum:
                    self.plotted_hist_momentum.Delete()
                    self.plotted_hist_momentum = None
                if self.fit_vertex:
                    self.fit_vertex.Delete()
                    self.fit_vertex = None
                if self.fit_momentum:
                    self.fit_momentum.Delete()
                    self.fit_momentum = None

                if self.c_vertex:
                    self.c_vertex.Close()
                    self.c_vertex = None
                if self.c_momentum:
                    self.c_momentum.Close()
                    self.c_momentum = None

                ROOT.gSystem.ProcessEvents()
                break
            else:
                print("Invalid input. Please enter 0, 1, or 2.")

def main():
    output_file = "flat_particle_ascii.hepmc" # Edit name of output file here
    n_events = 10000 # Edit the number of events here

    simulator = EventSimulator(output_file, n_events)
    simulator.generate_events()

    plotter = PlotManager(output_file)
    plotter.interactive_plot()

if __name__ == "__main__":
    main()























