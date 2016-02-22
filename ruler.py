## 2015-09-29

import openbabel, logging

global obconversion
obconversion = openbabel.OBConversion() 

global log
#logging.basicConfig(level=logging.WARNING)
logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)
log = logging

def replace_file_ext(matrix, new_ext, fieldname_prefix=""):
  for pos in matrix:
      matrix[pos][fieldname_prefix+"filename"]=matrix[pos][fieldname_prefix+"filename"][:matrix[pos][fieldname_prefix+"filename"].rfind(".")]+"."+new_ext

def read_files_create_first_matrix(dir_to_be_listed, file_startswith, file_endswith, fieldname_filename="filename"):
  import os
  matrix={}
  for file in os.listdir(dir_to_be_listed):
    if file.startswith(file_startswith) and file.endswith(file_endswith):
      filename=dir_to_be_listed+"/"+file
      index1=len(dir_to_be_listed+"/"+file_startswith)
      index2=filename.find(file_endswith,index1)
      pos = str(filename[index1:index2])
      bp=pos[:pos.find("_")]
      orientation=pos[pos.find("_")+1:]
      if orientation == "up":
	orientation=+1
      elif orientation=="down":
	orientation=-1
      log.debug(pos)
      matrix.update( { pos: {fieldname_filename: filename, "bp" : bp , "orientation" : orientation} } ) 
  return matrix

def sort_crossmatrix_into_bins(cm, number=10, start=0, stop=20, step=1):
  # if start is not set properly, glob_threshold optimization will not work... -> set to min(distance), if necessary...
  
  bins=range(start,stop,step)

  best={}
  for b in bins:
    best.update( { b : {"threshold" : None} } )

  glob_threshold = None
  
  for pos in cm:
    #print pos 
    # skip this entry if it is forbidden
    if cm[pos]["forbidden"] != 0:
      continue

    # get often-used values
    dist=cm[pos]["distance"]
    E = cm[pos]["sum_Etotal"]

    # if dist out of bin range: skip
    if dist > stop:
      continue
    
    # if all bins are full, check can be reduced to global threshold
    if glob_threshold != None:
      # if E bigger than global threshold: skip
      if E > glob_threshold:
	continue
    
    
    # identify corresponding bin
    for upperlimit in bins:
      if dist < upperlimit:
	# found corresponding bin
	
	# check if bin is not filled yet = no threshold.
	if best[upperlimit]["threshold"] == None:
	  # just add yourself to this bin
	  best[upperlimit].update({pos : cm[pos]})
	
	# obviously bin is already full.
	else:
	  # are you better than the threshold? than add yourself, no matter what.
	  if E <= best[upperlimit]["threshold"]:
	    best[upperlimit].update({pos : cm[pos]})
	  
	  # no, you're not. bye-bye!
	  else: 
	    break
	
	# was the bin filled right now?
	if len(best[upperlimit]) == number:
	  # this bin is full now; let's check if other bins are full as well...
	  print "bin "+str(upperlimit)+" was filled right now..." 
	  check = True
	  ## check if other bins are full, too:
	  for i in best:
	    if len(best[i]) < number:
	      #print "bin "+str(i)+" is not full yet..."
	      check = False
	      break
	  print check
	  # if so, set glob_threshold to maximum of thresholds
	  if check == True:
	    # set current threshold first...
	    best[upperlimit]["threshold"] = max([best[upperlimit][x]["sum_Etotal"] for x in best[upperlimit] if x != "threshold"])

	    # set global threshold to max of all thresholds
	    glob_threshold = max([best[i]["threshold"] for i in best])

	  
	# the bin is overfull! remove worst
	elif len(best[upperlimit]) > number:
	  old_threshold = best[upperlimit]["threshold"]
	  
	  for x in best[upperlimit]:
	    # catch threshold value
	    if x == "threshold":
	      continue
	    elif best[upperlimit][x]["sum_Etotal"] == best[upperlimit]["threshold"]:
	      del best[upperlimit][x]
	      break
	  
	  # the threshold can be set to the new worst...
	  best[upperlimit]["threshold"] = max([best[upperlimit][x]["sum_Etotal"] for x in best[upperlimit] if x != "threshold"])
	  
	  # was this maybe the global threshold? set it then, accordingly...
	  if glob_threshold != None and glob_threshold == old_threshold:
	    glob_threshold = max([best[i]["threshold"] for i in best])
	    

	# break this loop; continue with next entry
	break

  # clean up after me...
  for i in best:
    del(best[i]["threshold"])
	
  return best
  
def append_crossmatrix_with_shls(crossmatrix):
  SHL = {
    "-7" : [-75,-65],
    "-6" : [-65,-55],
    "-5" : [-55,-45],
    "-4" : [-45,-35],
    "-3" : [-35,-25],
    "-2" : [-25,-15],
    "-1" : [-15,-05],
    "0" : [-05,05],
    "+1" : [05,15],
    "+2" : [15,25],
    "+3" : [25,35],
    "+4" : [35,45],
    "+5" : [45,55],
    "+6" : [55,65],
    "+7" : [65,75],    
    }
    
  SHG = {
    "1" : SHL["-1"]+SHL["+7"],
    "2" : SHL["-2"]+SHL["+6"],
    "3" : SHL["-3"]+SHL["+5"],
    "4" : SHL["-4"]+SHL["+4"],
    "5" : SHL["-5"]+SHL["+3"],
    "6" : SHL["-6"]+SHL["+2"],
    "7" : SHL["-7"]+SHL["+1"],
    
    }
  
  for index in crossmatrix:
    #print index
    if index[0] > index[1]:
      current=index[0]
      compare=index[1]
    else:
      current=index[1]
      compare=index[0]
    #print current, compare
    bp_current=current[:current.find("_")]
    orientation_current=current[current.find("_")+1:]
    bp_compare=compare[:compare.find("_")]
    orientation_compare=compare[compare.find("_")+1:]

    crossmatrix[index].update( { "bp_current" : bp_current , "orientation_current" : orientation_current , "bp_compare" : bp_compare , "orientation_compare" : orientation_compare } )
    
    if float(crossmatrix[index]["bp_current"]) < float(crossmatrix[index]["bp_compare"]):
      small=float(crossmatrix[index]["bp_current"])
      big=float(crossmatrix[index]["bp_compare"])
    else:
      small=float(crossmatrix[index]["bp_compare"])
      big=float(crossmatrix[index]["bp_current"])
      
    for groove in SHG:
      if SHG[groove][0]<=small<SHG[groove][1] and SHG[groove][2]<=big<SHG[groove][3]:
	crossmatrix[index].update({ "SHG" : groove })
	break
    if "SHG" not in crossmatrix[index].keys():
      crossmatrix[index].update({ "SHG" : 0 })
	
    for locus in SHL:
      if SHL[locus][0]<=small<SHL[locus][1]:
	crossmatrix[index].update({ "SHL_small" : locus })
      if SHL[locus][0]<=big<SHL[locus][1]:
	crossmatrix[index].update({ "SHL_big" : locus })
    
    
  
def return_cross_matrix_for_binder_connections(matrix, binder_smiles, fieldname_prefix=""):
  storage=""
  ligands = []
  cross_matrix = {}
  for pos in matrix:
    filename = matrix[pos][fieldname_prefix+"filename"]
    ligands.append(openbabel.OBMol())
    obconversion.SetInFormat(filename[filename.rfind(".")+1:])
    if obconversion.ReadFile(ligands[-1], filename):
      log.info('Successfully read ligand '+str(filename))
    else:
      log.warning('Failed to read ligand '+str(filename))
    matrix[pos].update( {fieldname_prefix+"binder_atom" : get_binder_atom(ligands[-1], binder_smiles) } )
    #print pos
  used_indeces=[]
  for current in matrix:
    #print used_indeces
    used_indeces.append(current)
    log.debug("current "+current)
  
    if matrix[current][fieldname_prefix+"binder_atom"] != None:
      for compare in matrix:
        log.debug("compare "+compare)
        if (matrix[compare][fieldname_prefix+"binder_atom"] != None) and (compare not in used_indeces): 
          log.debug("new value added...")
          

          
          cross_matrix.update({(current, compare) : { "distance" : matrix[current][fieldname_prefix+"binder_atom"].GetDistance(matrix[compare][fieldname_prefix+"binder_atom"]), "sum_Etotal" : float(matrix[current][fieldname_prefix+"E_total"])+float(matrix[compare][fieldname_prefix+"E_total"]) } } )
          
          orientation = [current.split("_")[1], compare.split("_")[1]]
          number = [float(current.split("_")[0]), float(compare.split("_")[0])]
          
          length_binder_same_orientation = 4
          length_binder_opposite_orientation = 5
          
          if (orientation[0] == orientation[1]) and (abs(number[0]-number[1]) < length_binder_same_orientation):
	    cross_matrix[(current, compare)].update( { "forbidden" : 1 } ) #most likely overlapping
	  elif (orientation[0] != orientation[1]) and (abs(number[0]-number[1]) < length_binder_opposite_orientation):
	    cross_matrix[(current, compare)].update( { "forbidden" : 2 } ) #most likely overlapping
	  elif (orientation[0] != orientation[1]) and (abs(number[0]-number[1]) == length_binder_opposite_orientation):
	    cross_matrix[(current, compare)].update( { "forbidden" : 0.5 } ) #most likely neighboring
	  else:
	    cross_matrix[(current, compare)].update( { "forbidden" : 0 } ) #most likely okay
    else:
      log.info(str(current)+" did not contain a binder atom...")
  del ligands
  return cross_matrix
  
def append_matrix_with_center_of_geometry(matrix, fieldname_cog="cog", fieldname_prefix=""):
  import numpy
  ligands= []
  for pos in matrix:
    filename = matrix[pos][fieldname_prefix+"filename"]
    ligands.append(openbabel.OBMol())
    obconversion.SetInFormat(filename[filename.rfind(".")+1:])
    if obconversion.ReadFile(ligands[-1], filename):
      log.info('Successfully read ligand '+str(filename))
    else:
      log.warning('Failed to read ligand '+str(filename))
        
    cog = openbabel.vector3()
    for atom in openbabel.OBMolAtomIter(ligands[-1]):
      W = atom.GetVector()
      cog += W
    cog /= ligands[-1].NumAtoms()
    
    matrix[pos].update({ fieldname_prefix+fieldname_cog : cog })
  return
    
def append_matrix_with_delta_cog_to_previous_element(matrix, iterator, fieldname_deltacog = "delta_COG_to_previous"):
  import numpy
  foobar=None
  prev_cog=None
  for pos in iterator:
    if foobar==None:
      prev_cog=matrix[pos]["cog"]
      foobar=1 # skip first
      continue
    else:
      matrix[pos].update( { fieldname_deltacog : numpy.sqrt(matrix[pos]["cog"].distSq(prev_cog)) } )
      prev_cog=matrix[pos]["cog"]
  return
  
def append_matrix_with_delta_cog_in_same_line(matrix, field1="", field2="", fieldname_deltacog = "delta_COG_to_previous"):
  import numpy
  for pos in matrix:
    matrix[pos].update( { fieldname_deltacog : numpy.sqrt(matrix[pos][field2].distSq(matrix[pos][field1])) } )
  return
  
def append_matrix_with_content_of_file(matrix, fieldname_content = "RMSD", fieldname_prefix=""):
  for pos in matrix:
    filename = matrix[pos]["filename"]
    with open(filename, "r") as f:
      matrix[pos].update( { fieldname_prefix+fieldname_content : '"' + f.read() + '"'})
  
def append_matrix_with_rmsd_to_previous_element(matrix, iterator, fieldname_rmsd = "RMSD_to_previous"):
  import numpy
  ligands=[]
  foobar=None
  rmsd=0
  
  for pos in iterator:
    filename = matrix[pos]["filename"]
    ligands.append(openbabel.OBMol())
    obconversion.SetInFormat(filename[filename.rfind(".")+1:])
    if obconversion.ReadFile(ligands[-1], filename):
      log.info('Successfully read ligand '+str(filename))
    else:
      log.warning('Failed to read ligand '+str(filename))
    
    if foobar==None:
      foobar=1 # skip first
      continue
    else:
      #print "hi"
      foobar = openbabel.OBMolAtomIter(ligands[-2])
      for atom in openbabel.OBMolAtomIter(ligands[-1]):
	V = foobar.__next__().GetVector()
	W = atom.GetVector()
	rmsd+=(V.distSq(W))
      rmsd/=ligands[-1].NumAtoms()
      rmsd=numpy.sqrt(rmsd)
      
      matrix[pos].update( { fieldname_rmsd: rmsd } )
  return
    
  
def append_matrix_with_rmsd_in_same_line(matrix, prefix1="", prefix2="", fieldname_rmsd = "RMSD_to_previous"):
  tmp={}
  for pos in matrix:
    tmp.update( { pos: { prefix1+"filename" : matrix[pos][prefix1+"filename"], prefix2+"filename": matrix[pos][prefix2+"filename"] } } ) 
  with open(".pythonscript_to_chimera.py", "w") as f:
    f.write("matrix="+repr(tmp)+"\n")
    f.write("prefix1="+repr(prefix1)+"\n")
    f.write("prefix2="+repr(prefix2)+"\n")
    f.write("fieldname_rmsd="+repr(fieldname_rmsd)+"\n")
    f.write("""
import Midas
import chimera
  
for pos in matrix:
  filename1 = matrix[pos][prefix1+"filename"]
  filename2 = matrix[pos][prefix2+"filename"]

  Midas.open(filename1)
  Midas.open(filename2)
  matrix[pos].update( { "rmsd": Midas.rmsd("#0", "#1") } )
  chimera.runCommand("close #0 #1")

with open("/tmp/.pythonoutput_from_chimera.tmp", "w") as f:
  f.write("tmp2="+repr(matrix))
f.closed
  """)
  
  import os
  os.system("chimera --nogui --script .pythonscript_to_chimera.py")

  reimport={}
  execfile("/tmp/.pythonoutput_from_chimera.tmp", reimport)

  #print reimport["tmp2"]["+45.0_up"]["rmsd"]
  for pos in reimport["tmp2"]:
    matrix[pos].update( { fieldname_rmsd : reimport["tmp2"][pos]["rmsd"] } )
    
  #print matrix
  return
    
    
def append_matrix_with_linker_distances(matrix, linker_smiles):
  storage=""
  ligands = []
  for pos in matrix:
    filename = matrix[pos]["filename"]
    ligands.append(openbabel.OBMol())
    obconversion.SetInFormat(filename[filename.rfind(".")+1:])
    if obconversion.ReadFile(ligands[-1], filename):
      log.info('Successfully read ligand '+str(filename))
    else:
      log.warning('Failed to read ligand '+str(filename))
    linker_atoms=get_linker_atoms(ligands[-1],linker_smiles,linker_smiles)
    if len(linker_atoms)!=2:
      log.warning(file+" does not contain two recognizable linker ends!")
      storage+=filename+"\n"
    else:
      if (linker_atoms[0] != None) and (linker_atoms[1] != None): 
	matrix[pos].update({ "distance" : linker_atoms[0].GetDistance(linker_atoms[1]) } )
      else:
	log.warning("Uh-uh. Check this: "+str(matrix[pos]))
    if storage != "":
      log.warning("WARNING: The following files were corrupted / unrecognized:")
      log.warning(storage)
      log.warning("END OF WARNING")
  del ligands

def append_matrix_with_energies(matrix, fieldname_prefix=""):
  #import os
  
  how_to_recognize_energies={
  "E_total" : [ "Total", 2],
  "VdW" : ["__VdW", 1],
  "Lig-Prot" : ["Lig-Prot", 3],
  "Intramol" : ["Intramol", 4],
  }
  
  for energyform in how_to_recognize_energies:
    storage=""
    for pos in matrix:
      filename = matrix[pos][fieldname_prefix+"filename"]
      with open(filename, "r") as f:
	for line in f:
	  if how_to_recognize_energies[energyform][0] in line:
	    total_energy=line.split()[how_to_recognize_energies[energyform][1]]
	    break
      #print total_energy
      matrix[pos].update({ fieldname_prefix+energyform : total_energy})
      if storage != "":
	log.warning("WARNING: The following files were corrupted / unrecognized:")
	log.warning(storage)
	log.warning("END OF WARNING")

def write_matrix(matrix, output_filename, csv_delimiter = ","):
  csv = "index"
  if len(matrix.keys()) == 0:
    log.warning("Matrix for "+output_filename+" was empty. Did not write anything.") 
    return
  for index in matrix[matrix.keys()[-1]]:
    csv += csv_delimiter + str(index)
  csv += "\n"
    
  # 3.2. writes table lines  
  for row in sorted(matrix):
    # add matrix index (=pos)
    if type(row)==list or type(row)==tuple:
      csv += "("+" & ".join(row)+")" ##is this clever? ...
    else:
      csv += str(row)
    if type(matrix[row]) == dict:
      for column in matrix[row]:
        csv += csv_delimiter + str(matrix[row][column])# + csv_delimiter
    else:
      csv += csv_delimiter + str(matrix[row]) #+ csv_delimiter
      log.warning("Bad matrix format!") 
    csv += "\n"
    
  log.debug(csv)
  
  f = open(output_filename, 'w')
  f.write(csv)
  f.close()
  
  log.info("I have written the matrix to "+str(output_filename)+", right as you asked me to do.")



def get_linker_atoms(linker, smart1, smart2):
  ## workaround for symmetric linker atoms
  foobar = 0
  a1, a2 = None, None
  i=1
  for atom in openbabel.OBMolAtomIter(linker):
    log.debug(i)
    i+=1
    if foobar == 0 and atom.MatchesSMARTS(smart1):
      a1 = atom
      log.debug("found a1!")
    if atom.MatchesSMARTS(smart2) and atom != a1:
      a2 = atom
      log.debug("found a2!")
    elif atom.MatchesSMARTS(smart2) and atom == a1:
      foobar = 1
  return a1, a2

def get_binder_atom(binder, smarts):
  log.debug("Counting atoms to find the binder atom...")
  i=1
  for atom in openbabel.OBMolAtomIter(binder):
    log.debug(i)
    i+=1
    if atom.MatchesSMARTS(smarts):
      log.debug("found ya!")
      return atom

def connect_two_binders(binder1, binder1_smiles, binder2, binder2_smiles, linker, linker_smiles1, linker_smiles2):
  b1 = get_binder_atom(binder1, binder1_smiles).GetIdx()
  b2 = get_binder_atom(binder2, binder2_smiles).GetIdx()
  linker_atom1, linker_atom2 = get_linker_atoms(linker, linker_smiles1, linker_smiles2)
  l1 = linker_atom1.GetIdx()
  l2 = linker_atom2.GetIdx()

  linked = openbabel.OBMol()
  linked += binder1
  linked += linker
  linked += binder2
  
  idx_b1 = b1
  idx_l1 = binder1.NumAtoms() + l1
  idx_l2 = binder1.NumAtoms() + l2
  idx_b2 = binder1.NumAtoms() + linker.NumAtoms() + b2
  
  
  bond_order = 1
  linked.AddBond(idx_b1, idx_l1, bond_order)
  linked.AddBond(idx_b2, idx_l2, bond_order)

  foo = []
  for index in idx_b1, idx_l1, idx_l2, idx_b2:
    atom = linked.GetAtom(index)
    for child in openbabel.OBAtomAtomIter(atom):
      if child.IsHydrogen(): 
        j = child
    foo.append(j)
  for killme in foo:
    linked.DeleteAtom(killme)

  return linked
