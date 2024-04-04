#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'RunDD4hepInSerial.rb'
# Derek Anderson
# 03.21.2024
#
# Runs a series of DD4hep (ddsim/npsim) jobs. Two options for running:
#   type = :gun   -- run particle gun with parameters specified in
#                    provided steering file
#   type = :hepmc -- run specified hepmc file through DD4hep
#
# When providing command-line options (either via condor or locally), the
# order of arguments is:
#
#   ./RunDD4hepInSerial.rb
#     <num-iterations> <run-label> <steering-file> <input-file>
# -----------------------------------------------------------------------------

require 'fileutils'



# main body of script ---------------------------------------------------------

END {

  # default and configurable parameters ---------------------------------------

  # i/o parameters
  type    = :gun
  numiter = ARGV[0]
  input   = "forJetInitFix.e10h11pipXpim.d14m3y2024.hepmc"
  output  = "testOutOnLocal.d22m3y2024.edm4hep.root"
  steerer = ARGV[2]

  # input, output, and runnning directories
  in_dir  = "/sphenix/user/danderson/"
  out_dir = "/sphenix/user/danderson/"
  run_dir = "/sphenix/u/danderson/eic"

  # i/o parameters for condor running
  out_prefix = "testOutOnCondor"
  out_label  = ARGV[1]
  out_index  = "file"
  out_suffix = ".d22m3y2024.edm4hep.root"

  # fixed parameters ----------------------------------------------------------

  # job parameters
  numevts = 1000
  compact = "/opt/detector/epic-nightly/share/epic/epic.xml"

  # environment parameters
  exec  = "npsim"
  shell = "/sphenix/u/danderson/scripts/eic-shell"
  setup = "/opt/detector/setup.sh"

  # script names
  runner    = "RunDD4hepInShell.sh"
  simulator = "DoDD4hepSimulation.sh"

  # parse and run -------------------------------------------------------------

  iter = 0
  while iter < numiter.to_i

    # parse options
    parser = ArgParser.new
    parser.configure(
      args:  ARGV,
      steer: steerer,
      label: out_label,
      pref:  out_prefix,
      suff:  out_suffix,
      name:  output,
      sim:   simulator,
      run:   runner,
      index: out_index,
      count: iter
    )
    parser.parse()

    # run simulation
    handler = SimHandler.new
    handler.configure(
      outname:  parser.out_name,
      outdir:   out_dir,
      intype:   type,
      infile:   input,
      indir:    in_dir,
      optnum:   numevts,
      optdet:   compact,
      optsteer: parser.out_steer,
      sysexec:  exec,
      syssim:   parser.out_sim,
      sysrun:   parser.out_run,
      sysshell: shell,
      sysinit:  setup,
      rundir:   run_dir
    )
    handler.run()

    # increment counter
    iter = iter + 1

  end  # end while loop

}  # end main body of script



# class to parse input and return sensible job parameters ---------------------

class ArgParser

  # input
  attr_accessor :in_args
  attr_accessor :in_steer
  attr_accessor :in_label
  attr_accessor :in_pref
  attr_accessor :in_suff
  attr_accessor :in_name
  attr_accessor :in_sim
  attr_accessor :in_run
  attr_accessor :in_index
  attr_accessor :in_count

  # output
  attr_reader :out_steer
  attr_reader :out_label
  attr_reader :out_name
  attr_reader :out_sim
  attr_reader :out_run



  # external methods ----------------------------------------------------------

  def configure(
    args:  [],
    steer: "",
    label: "",
    pref:  "",
    suff:  "",
    name:  "",
    sim:   "DoDDSim.sh",
    run:   "RunShell.sh",
    index: out_num,
    count: iter
  )

    @in_args  = args
    @in_steer = steer
    @in_label = label
    @in_pref  = pref
    @in_suff  = suff
    @in_name  = name
    @in_sim   = sim
    @in_run   = run
    @in_index = index
    @in_count = count
    return

  end  # end :configure



  def parse()

    @out_steer = pick_steering_file()
    @out_label = pick_label()
    @out_name  = construct_output_name(@out_label)
    @out_sim   = construct_script_name(@in_sim, @out_label)
    @out_run   = construct_script_name(@in_run, @out_label)
    return

  end  # end :parse



  # internal methods ----------------------------------------------------------

  private

    def pick_steering_file()

      # use steering file from arguments if available
      steer = ""
      if @in_args.empty?
        steer = @in_steer.clone
      else
        steer = @in_args[2].clone
      end
      return steer

    end  # end :pick_steering_file



    def pick_label()

      # use label from arguments if available
      label = ""
      if @in_args.empty?
        label = @in_label.clone
      else
        label = @in_args[1].clone
      end

      # if iteration # provided, add to label
      if not @in_index.nil?
        label = label + @in_index + @in_count.to_s
      end
      return label

    end  # end :pick_label



    def construct_output_name(label)

      # first determine if output name is or should be empty
      name = ""
      if @in_args.empty?
        name = @in_name.clone
      end

      # then if output name is empty, construct one
      #   or append label if index name provided 
      if name.empty?
        name = @in_pref + "." + label + "." + @in_suff
        name.gsub!("..", ".")
      elsif not @in_index.nil?
        ending = "." + label + ".edm4hep.root"
        name.gsub!(".edm4hep.root", ending)
        name.gsub!("..", ".")
      end
      return name

    end  # end :construct_output_name



    def construct_script_name(input, label)

      # check if provided name ends in ".sh"
      script = input.clone
      if script.end_with?(".sh") and not label.empty?
        script.gsub!(".sh", ".#{label}.sh")
      end
      return script

   end  # end "construct_script_name

end  # end ArgParser



# class to consolidate options and run simulation -----------------------------

class SimHandler

  # getters
  attr_reader :out_file
  attr_reader :out_dir
  attr_reader :out_path
  attr_reader :in_type
  attr_reader :in_file
  attr_reader :in_dir
  attr_reader :in_orig
  attr_reader :in_path
  attr_reader :opt_nevt
  attr_reader :opt_det
  attr_reader :opt_steer
  attr_reader :sys_exec
  attr_reader :sys_sim
  attr_reader :sys_run
  attr_reader :sys_shell
  attr_reader :sys_init
  attr_reader :sim_path
  attr_reader :run_path
  attr_reader :run_dir



  # external methods ----------------------------------------------------------

  public

    def configure(
      outname:  "",
      outdir:   "./",
      intype:   :gun,
      infile:   nil,
      indir:    "./",
      optnum:   1,
      optdet:   "epic.xml",
      optsteer: "./steer.py",
      sysexec:  "ddsim",
      syssim:   "DoSim.sh",
      sysrun:   "RunShell.sh",
      sysshell: "~/scripts/eic-shell",
      sysinit:  "~/scripts/initialize-eic-detectors",
      rundir:   "./"
    )

      @out_file  = outname
      @out_dir   = outdir
      @in_type   = intype
      @in_file   = infile
      @in_dir    = indir
      @opt_nevt  = optnum
      @opt_det   = optdet
      @opt_steer = optsteer
      @sys_exec  = sysexec
      @sys_sim   = syssim
      @sys_run   = sysrun
      @sys_shell = sysshell
      @sys_init  = sysinit
      @run_dir   = rundir
      return

    end  # end :configure



    def run()

      # parse what was given to handler
      prepare_for_running()

      # make job scripts to run
      create_sim_script(@sim_path) if not File.exists?(@sim_path)
      create_run_script(@run_path) if not File.exists?(@run_path)

      # run simulation
      system("#{@sys_shell} -- #{@run_path}")

      # remove job scripts
      clean_up()
      return

    end  # end :run



  # internal methods ----------------------------------------------------------
  private

    def prepare_for_running()

      # construct paths to sim/run scripts and output
      @sim_path = @run_dir + "/" + @sys_sim
      @run_path = @run_dir + "/" + @sys_run
      @out_path = @run_dir + "/" + @out_file
      @sim_path.gsub!("//", "/")
      @run_path.gsub!("//", "/")
      @out_path.gsub!("//", "/")
      @sim_path.gsub!("..", ".")
      @run_path.gsub!("..", ".")
      @out_path.gsub!("..", ".")

      # if needed make sure input is in run directory
      if in_type == :hepmc

        # construct paths
        @in_orig = @in_dir  + "/" + @in_file
        @in_path = @run_dir + "/" + @in_file
        @in_orig.gsub!("..", ".")
        @in_path.gsub!("//", "/")
        @in_orig.gsub!("//", "/")
        @in_path.gsub!("..", ".")

        # copy to run dir
        if not File.exists?(@in_path)
          FileUtils.cp(@in_orig, @in_path)
        end
      end
      return

    end  # end :prepare_for_running



    def create_sim_script(path)

      # script input
      args_array = ["setup", "rundir", "steerer", "compact", "numevts", "output"]
      args_array << "input" if in_type == :hepmc

      # create block to assign arguments
      args_block = ""
      args_array.each_with_index do |arg, iArg|
        args_num    = iArg + 1
        args_block += arg
        args_block += "=$"
        args_block += "#{args_num}"
        args_block += "\n"
      end

      # create block for setup
      setup_block  = ["source $setup", "cd $rundir"].join("\n")
      setup_block += "\n"

      # create command
      command = ""
      if in_type == :gun
        command += "#{@sys_exec} --steeringFile $steerer --compactFile $compact -G -N $numevts --outputFile $output"
      elsif in_type == :hepmc
        command += "#{@sys_exec} --steeringFile $steerer --compactFile $compact -I $input -N $numevts --outputFile $output"
      end
      command += "\n"

      # create script
      body = [args_block, setup_block, command].join("\n")
      File.write(path, body)
      FileUtils.chmod("u+x", path)
      return

    end # end :create_sim_script



    def create_run_script(path)

      # arguments to give to script
      args_array = [@sys_init, @run_dir, @opt_steer, @opt_det, @opt_nevt, @out_path]
      args_array << @in_path if in_type == :hepmc

      # create command to run
      arguments = args_array.join(" ")
      command   = ["source", @sim_path, arguments].join(" ")

      File.write(path, command)
      FileUtils.chmod("u+x", path)
      return

    end  # end :create_run_script



    def clean_up()

      # move output to output directory
      FileUtils.mv(@out_path, @out_dir)

      # remove sim/run scripts
      FileUtils.rm @sim_path
      FileUtils.rm @run_path

      # delete input if needed
      FileUtils.rm @in_path if in_type == :hepmc
      return

    end  # end :clean_up

end  # end SimHandler

# end -------------------------------------------------------------------------
