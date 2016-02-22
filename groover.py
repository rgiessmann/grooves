#!/usr/bin/env python

# version 2015-09-16
## needs DEBUGging for: 
## -placing ligands to / around bp0 -> only partial recognition of bp0 
## -generating dummies (from where to where, when? helical / minorgroove)

import openbabel
import logging
import sys
import getopt

### +++ Start of Default Settings +++ ###

ncp_format = "pdb"

ligand_format = "pdb"


dna_acceptable_nucleobase_residues = ['DA', 'DC', 'DG', 'DT', "A", "C", "G", "T"]
logging.basicConfig(level=logging.WARNING)
allow_interactive_mode = False

### +++ End of Settings +++ ###
###  don't edit below here! ###

global log
log = logging

### +++ Get User Settings from Commandline +++ ###

#-h, --help:

#-i, --iterate [False]
#    --iterate-from=
#    --iterate-to=
#    --iter-step= [1]

#-a, --align-to= [False]
#    --down-pose [False]
#    --both-poses [False]

#-g, --generate-all-dummies [True]
#    --generate-dummies-from=
#    --generate-dummies-to=
#-f, --finegrain-step= [1]
#-k, --keep-dummies [False]

#    --dna-leading-strand= [I]
#    --dna-lagging-strand= [J]
#    --dna-leading-start= [-72]
#    --dna-lagging-start= [+72]

#    --bp-counting-parallel [False]
#    --minor-groove-offset= [2]

#-o, --output-prefix= [output/]

#[1]   => dna-file
#[2]   => ligand-file

# ++ Start of Program ++ #

def main(argv=""):
  # Set default variables
  iterate = False
  iterate_from = None
  iterate_to = None
  iter_step = 1
  
  align_to = None
  
  up_pose = True
  both_poses = False
  
  generate_all_dummies = True
  generate_dummies_from = None
  generate_dummies_to = None

  global finegrain_step
  finegrain_step = 1
  
  global remove_dummies
  remove_dummies = True
  
  global dna_leading_strand_basepairs_in_advance_of_lagging_strand_when_looking_at_minor_groove
  dna_leading_strand_basepairs_in_advance_of_lagging_strand_when_looking_at_minor_groove = 2
  
  global dna_bp_counting_antiparallel
  dna_bp_counting_antiparallel = 1

  global dna_leading_strand
  dna_leading_strand = 'I'
  global  dna_leading_strand_starts_at_position
  dna_leading_strand_starts_at_position = -72
  global  dna_lagging_strand
  dna_lagging_strand = 'J'
  global  dna_lagging_strand_position_at_leading_strand_start
  dna_lagging_strand_position_at_leading_strand_start = +72


  output_prefix = "output/"
  
  # check given command-line arguments
  try:
    opts, remaining_args = getopt.getopt(argv,"dhia:gf:o:k",["help","iterate","iterate-from=","iterate-to=","iter-step=","align-to=","down-pose","both-poses","generate-all-dummies","generate-dummies-from=","generate-dummies-to=","finegrain-step=","keep-dummies","dna-leading-strand=","dna-lagging-strand=","dna-leading-start=","dna-lagging-start=","bp-counting-parallel","minor-groove-offset=","output-prefix="])
  except getopt.GetoptError:
    print 'You provided unusual arguments. Call me with -h to learn more.'
    sys.exit(2)
  for opt, arg in opts:
    if opt in ('-h', '--help'):
      print 'The following options are available:'
      print """
      -h, --help:

      -i, --iterate [False]
         --iterate-from=
         --iterate-to=
         --iter-step= [1]

      -a, --align-to= [False]
         --down-pose [False]
         --both-poses [False]

      -g, --generate-all-dummies [True]
         --generate-dummies-from=
         --generate-dummies-to=
      -f, --finegrain-step= [1]
      -k, --keep-dummies [False]

         --dna-leading-strand= [I]
         --dna-lagging-strand= [J]
         --dna-leading-start= [-72]
         --dna-lagging-start= [+72]

         --bp-counting-parallel [False]
         --minor-groove-offset= [2]

      -o, --output-prefix= [output/]

      [1]   => dna-file
      [2]   => ligand-file"""
      sys.exit()
    elif opt in ("-d"):
      logging.basicConfig(level=logging.DEBUG)
      global log
      log = logging
    elif opt in ("-i", "--iterate"):
      #log.info("Using iteration mode for aligning of ligands...")
      iterate = True
    elif opt in ("--iterate-from"):
      iterate_from = float(arg)
      #log.info("Set iterate-from: "+iterate_from)
    elif opt in ("--iterate-to"):
      iterate_to= float(arg)
      #log.info("Set iterate-to: "+iterate_to)
    elif opt in ("--iter-step"):
      iter_step = float(arg)
    elif opt in ("-a", "--align-to"):
      center_ligand_at = float(arg)
      align_to = True
    elif opt in ("--down-pose"):
      up_pose = False
    elif opt in ("--both-poses"):
      both_poses = True
    elif opt in ('-g', '--generate-all-dummies'):
      generate_all_dummies = True
    elif opt in ("--generate-dummies-from"):
      generate_dummies_from = int(arg)
      generate_all_dummies = False
    elif opt in ("--generate-dummies-to"):
      generate_dummies_to = int(arg)
      generate_all_dummies = False
    elif opt in ("-f", "--finegrain-step"):
      arg = float(arg)
      if arg > 1:
	log.warning("You provided the finegrain-step argument with a value greater than 1. This is non-sense, and will be ignored => finegrain-step = 1.")
      elif int(1/arg) != 1/arg: #checks if finegrain-step is sensible, i.e. smooth divisor of 1
	log.warning("You provided the finegrain-step argument with a value which will not add up to 1. This will result in non-uniform distribution of interpolated dummies. I will continue anyway, interpolating with the given steps until reaching 1.")
	finegrain_step = float(arg)
      else:
	log.info("Using finegraining, i.e. interpolation between dummy atoms...")
	finegrain_step = float(arg)
    elif opt in ("--dna-leading-strand"):
        dna_leading_strand = str(arg)
    elif opt in ("--dna-lagging-strand"):
        dna_lagging_strand = str(arg)
    elif opt in ("--dna-leading-start"):
        dna_leading_strand_starts_at_position = int(arg)
    elif opt in ("--dna-lagging-start"):
        dna_lagging_strand_position_at_leading_strand_start = int(arg)
    elif opt in ("--bp-counting-parallel"):
        dna_bp_counting_antiparallel = False
    elif opt in ("--minor-groove-offset"):
      dna_leading_strand_basepairs_in_advance_of_lagging_strand_when_looking_at_minor_groove = int(arg)
    elif opt in ("-o", "--output-prefix"):
      output_prefix = str(arg)+"/"
    elif opt in ("-k", "--keep-dummies"):
      remove_dummies = False
  
  # check mandatory file fields
  if not len(remaining_args) == 2:
    log.critical("You did not provide input files. Call me with -h to learn more.")
    sys.exit()
  else:
    global dna_file
    dna_file = remaining_args[0]
    ligand_file = remaining_args[1]
      
  # check argument combinations; set advanced (=clever) default values
  if (iterate_from and not iterate_to) or (iterate_to and not iterate_from):
    log.critical("You provided only iterate-from or iterate-to, the other is missing. This doesn't make sense; quitting.")
    sys.exit(2)
  elif iterate and not iterate_to and not iterate_from:
    log.critical("You wanted to iterate, but did neither provide iterate-from nor iterate-to. Please do this; I am quitting.")
    sys.exit(2) 
  elif iterate_from and iterate_to:
    if iterate_from > iterate_to:
      log.info('You gave iterate-from and iterate-to in reverse order. I changed that for you...')
      _tmp = iterate_from
      iterate_from = iterate_to
      iterate_to = _tmp
    if not iterate:
      log.info('You forgot to set iterate. Setting iterate to True, as you provided iterate-from and iterate-to.')
      iterate = True
    if not iter_step:
      log.info('Setting iter_step to 1 [default].')
      iter_step = 1
  if not iterate and not align_to:
    log.info("You chose not to perform any ligand alignments. Alright, then...")
    no_alignment = True
  elif align_to and iterate:
    log.info("You wanted to do an iteration and a single alignment. Don't know what to do; quitting.")
    sys.exit(2)
  else:
    no_alignment = None
     
     
  # ... done with configurations; start the real work!
  
  ## create ncp object
  ncp = openbabel.OBMol()
  ## create ligand object
  ligand = openbabel.OBMol()
  ## read files
  ## WATCH OUT: works only for single molecule files!
  obconversion = openbabel.OBConversion()
  obconversion.SetInFormat(ncp_format)
  if obconversion.ReadFile(ncp, dna_file):
    log.info('Successfully read DNA containing file ' + str(dna_file))
  obconversion.SetInFormat(ligand_format)
  if obconversion.ReadFile(ligand, ligand_file):
      log.info('Successfully read ligand file ' + str(ligand_file))
  
  leading_strand_phosphates, lagging_strand_phosphates = get_all_phosphates(ncp)
  
  if generate_all_dummies:
    # overriding even otherwisely set values...
    generate_dummies_from = min([int(i) for i in leading_strand_phosphates])+dna_leading_strand_basepairs_in_advance_of_lagging_strand_when_looking_at_minor_groove #-1 this does not work for antiparallel counting
    generate_dummies_to = max([int(i) for i in leading_strand_phosphates])-1
    log.debug("Generating dummies from "+ str(generate_dummies_from) + " to "+ str(generate_dummies_to))
  ncp_dummy_atoms, helix_dummy_atoms = generate_dummies(leading_strand_phosphates, lagging_strand_phosphates, generate_dummies_from, generate_dummies_to, finegrain_step)
  
  ## write out mol2 file
  obconversion.SetOutFormat("mol2")
  obconversion.WriteFile(ncp_dummy_atoms, output_prefix+"dummies.mol2")
  obconversion.SetOutFormat("pdb")
  obconversion.WriteFile(ncp_dummy_atoms, output_prefix+"dummies.pdb")
  obconversion.SetOutFormat("mol2")
  obconversion.WriteFile(helix_dummy_atoms, output_prefix+"helix.mol2")
  obconversion.SetOutFormat("pdb")
  obconversion.WriteFile(helix_dummy_atoms, output_prefix+"helix.pdb")
    
  # break here, if no alignment is wished
  if no_alignment:
    return
  
  if up_pose == True:
    posename = "up"
    antipose = "down"
  else:
    posename = "down"
    antipose = "up"
  
  if iterate:
    center = iterate_from
    iter_to = iterate_to
    iteration_step = iter_step
  if not iterate:
    center = center_ligand_at
    iteration_step = 1
    iter_to = center_ligand_at+1
    
  while center < iter_to: # + iteration_step):
    align_to_these_atoms = select_dummy_atoms(ncp_dummy_atoms, helix_dummy_atoms, center, up_pose)
    aligned_ligand, rmsd = align_ligand(align_to_these_atoms, ligand)
   
    if center > 0:
      sign = "+"
    else:
      sign = ""
    
    if aligned_ligand != None:
      log.debug("RMSD: "+ str(rmsd))
      obconversion.WriteFile(aligned_ligand, output_prefix+"ligand_aligned_to_bp" + sign + str(center) + "_"+posename+".pdb")
      with open(output_prefix+"ligand_aligned_to_bp" + sign + str(center) + "_"+posename+".log", "w") as f:
	f.write(str(rmsd)+"\n")
    
    if both_poses:
      align_to_these_atoms = select_dummy_atoms(ncp_dummy_atoms, helix_dummy_atoms, center, (not up_pose))
      aligned_ligand, rmsd = align_ligand(align_to_these_atoms, ligand)
      
      if aligned_ligand != None:
	obconversion.WriteFile(aligned_ligand, output_prefix+"ligand_aligned_to_bp" + sign + str(center) + "_"+antipose+".pdb")
	with open(output_prefix+"ligand_aligned_to_bp" + sign + str(center) + "_"+antipose+".log", "w") as f:
	  f.write(str(rmsd)+"\n")
	  
    center += iteration_step

  return

def select_dummy_atoms(all_dummies, helix_dummies, center, positive_pose=True):
  select = openbabel.OBMol()
  search_for = ""
  
  ## 1.1 start the search for the correct atoms
  log.info('Searching for dummies with center position ' + str(center) + "...")
  
  ## 1.2 find corresponding helix dummy
  search_for =  str(float(center)-1)
  log.debug("Starting search for: "+search_for)
  for res in openbabel.OBResidueIter(helix_dummies):
    log.debug(res.GetNumString())
    if float(res.GetNumString()) == float(search_for):
      log.info('... found dummy number 1 at helix.')
      select.AddResidue(res)
      for atom in openbabel.OBResidueAtomIter(res): 
	select.AddAtom(atom)
      break
  
  ## 1.2.1. find next helix dummy
  search_for =  str(float(center))
  for res in openbabel.OBResidueIter(helix_dummies):
    log.debug(res.GetNumString())
    if float(res.GetNumString()) == float(search_for):
      log.info('... found dummy number 2 at helix.')
      select.AddResidue(res)
      for atom in openbabel.OBResidueAtomIter(res): 
	select.AddAtom(atom)
      break
  
  ## 1.3 find next-to-center dna dummy
  search_for =  str(float(center)-1)
  for res in openbabel.OBResidueIter(all_dummies):
    log.debug(res.GetNumString())
    if float(res.GetNumString()) == float(search_for):
      log.info('... found dummy before center.')
      select.AddResidue(res)
      for atom in openbabel.OBResidueAtomIter(res): 
	select.AddAtom(atom)
      break
      
  ## 1.4 find center dna atom
  search_for =  str(float(center)+2)
  for res in openbabel.OBResidueIter(all_dummies):
    log.debug(res.GetNumString())
    if float(res.GetNumString()) == float(search_for):
      log.info('... found dummy after center.')
      select.AddResidue(res)
      for atom in openbabel.OBResidueAtomIter(res): 
	select.AddAtom(atom)
      break

  ## 1.3 check for success
  if select.NumAtoms() != 4:
    log.critical("Could not select DNA dummy atoms to align the ligand to. I have selected "+str(select.NumAtoms())+" atoms. Alignment will fail!")

  return select
  
  
def get_all_phosphates(ncp):
  # this function generates all dummy-atoms for the DNA minor groove

  ## 1. get positions of phosphates

  ### 1.1 generate phosphate storage variables as dictionaries
  dna_leading_strand_phosphates = {} #OLD: range(0,dna_leading_strand_length)
  dna_lagging_strand_phosphates = {} #OLD: range(0,dna_lagging_strand_length)

  ### 1.2 search all atoms for phosphates of the DNA strands
  for atom in openbabel.OBMolAtomIter(ncp):
    if atom.IsPhosphorus():
      # 1.2.1 we found a phosphate...
      _current_res=atom.GetResidue()
      if _current_res.GetChain() == dna_leading_strand and _current_res.GetName() in dna_acceptable_nucleobase_residues:
	# 1.2.1.1 we found a phosphate of a nucleobase in the leading chain; sort it to the right position (residue#) into the dict
	pos = _current_res.GetNumString().strip() #-dna_leading_strand_offset
	vec = atom.GetVector()
	log.debug('found phosphate in the leading chain at position '+str(pos)+'...')
	dna_leading_strand_phosphates.update({pos : vec})
      if _current_res.GetChain() == dna_lagging_strand and _current_res.GetName() in dna_acceptable_nucleobase_residues:
	# 1.2.1.2 we found a phosphate of a nucleobase in the lagging chain; sort it to the right position (residue#) into the list
	pos = _current_res.GetNumString().strip() #)-dna_leading_strand_offset
	vec = atom.GetVector()
	log.debug('found phosphate in the lagging chain at position '+str(pos)+'...')
	dna_lagging_strand_phosphates.update({pos : vec})
	
  #### 1.3 check phosphate recognition for completeness!
  #log.info('checking phosphates in leading strand for completeness...')
  #for _error_check_count in range(0,dna_leading_strand_length):
    #log.debug('... visiting position ' +str(_error_check_count) + " in leading strand")
    #if not _error_check_count in dna_leading_strand_phosphates: #OLD: .keys():
      #log.warning("Could not define all phosphates in the leading strand! Phosphate at position "+str(_error_check_count+dna_leading_strand_offset)+" is missing!")
  ##log.info('... done.')
  
  #log.info('checking phosphates in lagging strand for completeness...')
  #for _error_check_count in range(0,dna_lagging_strand_length):
    #log.debug('... visiting position ' +str(_error_check_count) + " in lagging strand")
    #if not _error_check_count in dna_lagging_strand_phosphates: #OLD: .keys():
      #log.warning("Could not define all phosphates in the lagging strand! Phosphate at position "+str(_error_check_count+dna_leading_strand_offset)+" is missing!")
  ##log.info('... done.')
  
  return dna_leading_strand_phosphates, dna_lagging_strand_phosphates

def generate_dummies(dna_leading_strand_phosphates, dna_lagging_strand_phosphates, start_at_position, stop_at_position, finegraining_step=1, *args):
  ## this function generates (all) dummy atoms from the given phosphates

  log.info('generating dummy atom coordinates from phosphates...')
  
  ## 0.9 create dummy atom object
  dummies = openbabel.OBMol()      
  dummies.SetTitle("minor groove dummy atom cloud for "+dna_file)
  
  helix = openbabel.OBMol()      
  helix.SetTitle("minor groove dummy atom cloud for "+dna_file)
  
  
  start = 0
  stop = 0
  step = 0
  
  # 1.0 check operation mode; set start and stop accordingly
  if start_at_position < stop_at_position: #construct_dummies_from_leading_strand_position > construct_dummies_to_leading_strand_position :
    log.info("Counting upwards for dummy generation...")
    step = +1
    start = start_at_position
    stop = stop_at_position+1
  elif start_at_position > stop_at_position:
    log.info("Counting downwards for dummy generation...")
    step = -1
    start = start_at_position
    stop = stop_at_position-1
    
  # 1.1 get coordinates of relevant phosphates
  log.info('Getting coordinates of phosphates ...')
  for i in range(start, stop, step):
    lead_index = str(i)
    log.debug(lead_index)
    progress = i - dna_leading_strand_starts_at_position
    log.debug(progress)
    for mode in ["minorgroove", "helix"]:
      
      if dna_bp_counting_antiparallel:
	if mode == "minorgroove":
	  lag_index = str(dna_lagging_strand_position_at_leading_strand_start - progress + dna_leading_strand_basepairs_in_advance_of_lagging_strand_when_looking_at_minor_groove)
	elif mode == "helix":
	  lag_index = str(dna_lagging_strand_position_at_leading_strand_start - progress)
      else:
	if mode == "minorgroove":
	  lag_index = str(dna_lagging_strand_position_at_leading_strand_start + progress - dna_leading_strand_basepairs_in_advance_of_lagging_strand_when_looking_at_minor_groove)
	elif mode == "helix":
	  lag_index = str(dna_lagging_strand_position_at_leading_strand_start + progress)
      
      log.debug(lag_index)  
      log.debug('... Visiting position ' + lead_index + ' in leading strand')
      
      _vector_leading = dna_leading_strand_phosphates[lead_index] 
      _vector_lagging = dna_lagging_strand_phosphates[lag_index]
      
      ## 1.1 create coordinates of dummy atoms with linear combination of vectors   
      _coord_dummy_atom = openbabel.vector3()
      _coord_dummy_atom.SetX((_vector_leading.GetX() + _vector_lagging.GetX())/2)
      _coord_dummy_atom.SetY((_vector_leading.GetY() + _vector_lagging.GetY())/2)
      _coord_dummy_atom.SetZ((_vector_leading.GetZ() + _vector_lagging.GetZ())/2)
      
      log.debug("... Coordinates of dummy: " + str(_coord_dummy_atom.GetX()) + " " + str(_coord_dummy_atom.GetY()) + " " + str(_coord_dummy_atom.GetZ()))
      
      ## 1.2 create atom representation in openbabel
      if mode == "minorgroove":
	_new_atom = dummies.NewAtom()
      elif mode == "helix":
	_new_atom = helix.NewAtom()
      _new_atom.SetType("Du")
      _new_atom.SetAtomicNum(0)
      _new_atom.SetVector(_coord_dummy_atom)
      
      ## 1.3 create and add to residue representation
      
      if mode == "minorgroove":
	_new_res = dummies.NewResidue()
	_new_res.SetName("DUM") 
      elif mode == "helix":
	_new_res = helix.NewResidue()
	_new_res.SetName("HDU") 
      _new_res.SetNum(str(float(i)))
      log.debug("Created atom #"+_new_res.GetNumString())
      
      _new_res.SetChain("W")
      #_new_res.SetChainNum(99)
      _new_res.AddAtom(_new_atom)
      _new_res.SetAtomID(_new_atom, "dummy")
      _new_res.SetHetAtom(_new_atom, 1)
      
      ## 1.4 interpolate between dummies, if wanted (= finegraining)
      if (finegraining_step < 1) and (i != stop-1): # construct_dummies_to_leading_strand_position ):
	mantissa = finegraining_step
		
	lead_index_next = str(int(lead_index)+step)
	log.debug(lead_index_next)
	if dna_bp_counting_antiparallel:
	  lag_index_next = str(int(lag_index)-step)
	else:
	  lag_index_next = str(int(lag_index)+step)
	log.debug(lag_index_next)
	_next_vector_leading = dna_leading_strand_phosphates[lead_index_next]
	_next_vector_lagging = dna_lagging_strand_phosphates[lag_index_next]
	while mantissa < 1-finegrain_step:
	  ## 1.4.1 create coordinates of dummy atoms with linear combination of vectors   
	  _coord_dummy_atom = openbabel.vector3()
	  _coord_dummy_atom.SetX((1-mantissa)*(_vector_leading.GetX() + _vector_lagging.GetX())/2+(mantissa)*(_next_vector_leading.GetX() + _next_vector_lagging.GetX())/2)
	  _coord_dummy_atom.SetY((1-mantissa)*(_vector_leading.GetY() + _vector_lagging.GetY())/2+(mantissa)*(_next_vector_leading.GetY() + _next_vector_lagging.GetY())/2)
	  _coord_dummy_atom.SetZ((1-mantissa)*(_vector_leading.GetZ() + _vector_lagging.GetZ())/2+(mantissa)*(_next_vector_leading.GetZ() + _next_vector_lagging.GetZ())/2)
	  log.debug("... Coordinates of dummy: " + str(_coord_dummy_atom.GetX()) + " " + str(_coord_dummy_atom.GetY()) + " " + str(_coord_dummy_atom.GetZ()))
	  
	  ## 1.4.2 create atom representation in openbabel
	  if mode == "minorgroove":
	    _new_atom = dummies.NewAtom()
	  elif mode == "helix":
	    _new_atom = helix.NewAtom()
	  _new_atom.SetType("Du")
	  _new_atom.SetAtomicNum(0)
	  _new_atom.SetVector(_coord_dummy_atom)
	  
	  ## 1.4.3 create and add to residue representation
	  if mode == "minorgroove":
	    _new_res = dummies.NewResidue()
	    _new_res.SetName("DUM") 
	  elif mode == "helix":
	    _new_res = helix.NewResidue()
	    _new_res.SetName("HDU") 
	  
	  if step > 0: 
	    _new_res.SetNum(str(i+mantissa))
	  else:
	    _new_res.SetNum(str(i-mantissa))
	  
	  log.debug("Created atom #"+_new_res.GetNumString())
	  
	  _new_res.SetChain("W")
	  #_new_res.SetChainNum(99)
	  _new_res.AddAtom(_new_atom)
	  _new_res.SetAtomID(_new_atom, "dummy")
	  _new_res.SetHetAtom(_new_atom, 1)
	  
	  ## 1.4.4 try if there is a next step to take...
	  mantissa += finegraining_step

  return dummies, helix    
  
  

def get_dummies(ligand):
  # reads ligand, extracts dummies
  log.info('Extracting dummy atoms from ligand ...')
  
  # 0.9 initialize output variable
  ligand_dummies = openbabel.OBMol()
  
  # 1.0 search for three dummy atoms in the ligand file
  for atom in openbabel.OBMolAtomIter(ligand):  
    #_current_res=atom.GetResidue()
    #print atom.GetAtomicNum()
    if atom.GetAtomicNum() == 0:
      ## 1.1 found a dummy atom by its atomic number
      ligand_dummies.AddAtom(atom)
      log.debug('... coordinates of extracted dummy atoms from ligand: ' + str(atom.GetVector().GetX()) + " " + str(atom.GetVector().GetY()) + " " + str(atom.GetVector().GetZ()))
  
  # 1.1 check for possibility to align
  if ligand_dummies.NumAtoms() != 4:
    log.critical("Could not find the expected number of dummy atoms in the ligand file, but it were "+str(ligand_dummies.NumAtoms())+" ... -- alignment is going to fail!")

  return ligand_dummies
  
  
def align_ligand(dummies, ligand):
  # fit dummy atoms of ligand to defined positions
  log.info('Aligning ligand dummy atoms to desired dummy atoms...')

  # 0.9 create local copy, as this function would otherwise modify the given ligand 
  aligned_ligand = openbabel.OBMol(ligand)
  
  # 1.0 get dummy atoms from ligand
  log.debug('... get dummy atoms of ligand')
  ligand_dummies = get_dummies(ligand)
  
  # 1.1 get translation vector from read-in position to origin
  log.info('... determing translation vector from read-in to origin')
  translation = ligand_dummies.Center(1)
  
  ## DEBUG
  #obconversion = openbabel.OBConversion()
  #obconversion.SetOutFormat("pdb")
  #obconversion.WriteFile(ligand_dummies,"ligand_dummies_centered.pdb")
  
  # 1.2 initialize OBAlign for alignment to final destination
  log.info('... doing the alignment for dummy atoms')
  aligner = openbabel.OBAlign(dummies, ligand_dummies)
  success=aligner.Align()
  
  if success == False:
    return None, None
    
  #log.info('... done.')
  rmsd=aligner.GetRMSD()
  log.debug('RMSD of alignment: ' + str(rmsd))

  ## 1.2.1 get Rotation Matrix for alignment 
  log.info('... determining the rotation matrix')
  rotation_matrix = aligner.GetRotMatrix()
  rot = openbabel.double_array([1,2,3,4,5,6,7,8,9])
  rotation_matrix.GetArray(rot)
  
  # .. only for debugging:
  arewedebugging = log.getLogger()
  if arewedebugging.isEnabledFor(logging.DEBUG):
    log.debug('--- rotation matrix ---')
    for i in range(0,9): log.debug(str(i)+ " : " + str(rot[i]))

  # 1.3 generate positioning vector 
  ## NB: centering would not work because of rotation!!
  ## update cooordinates to new value
  aligner.UpdateCoords(ligand_dummies)
  log.info('... generating positioning vector')
  positioning = openbabel.vector3()
  ## calculate the vector for positioning to destination
  n = 0
  for atom in openbabel.OBMolAtomIter(ligand_dummies):
    n += 1
    positioning += atom.GetVector()
  positioning /= n 

  # 1.4 move all ligand atoms to fit the pose the dummy atoms have been aligned to
  
  ## 1.4.1 generate inverted translation vector to translate ligand to origin
  translation_to_origin = translation
  translation_to_origin *= -1

  ## 1.4.2 generate inverted rotation matrix to reconstruct alignment results
  rotation_inversion = rotation_matrix.transpose()
  rot_inv = openbabel.double_array([1,2,3,4,5,6,7,8,9])
  rotation_inversion.GetArray(rot_inv)
    
  ## 1.4.3 apply translation to origin, rotation, and translation to final destination
  aligned_ligand.Translate(translation_to_origin)
  aligned_ligand.Rotate(rot_inv)
  aligned_ligand.Translate(positioning)
  
  ## 1.5 clean output ligand of dummy atoms, if desired
  if remove_dummies:
    log.info('Cleaning the ligand of unwanted dummies...')
    _temp_atom = []
    for atom in openbabel.OBMolAtomIter(aligned_ligand):
      #aligned_ligand.AddResidue(atom.GetResidue())
      #aligned_ligand.AddAtom(atom)
      if atom.GetAtomicNum() == 0:
	_temp_atom.append(atom)
    for a in _temp_atom:
      aligned_ligand.DeleteAtom(a)
  #else:
    #aligned_ligand = ligand

  log.info('... returning the aligned ligand.')
  return aligned_ligand, rmsd

  
  
## +++ MAIN SECTION: makes calls when script is started +++ ###
  
if allow_interactive_mode:
  log.info('Script started successfully.')
  log.info('Run program with "main()" if you like...')
elif __name__ == "__main__":
  main(sys.argv[1:])
