try:
    import pymol
    from pymol import cmd
    from pymol import stored
    from pymol import CmdException
except:
    cmd = None
    pymol = None
    print("Pymol is most likely not installed")

import os
import re
import itertools
from sys import prefix

# taken from http://stackoverflow.com/a/11301781/952600
try:
    basestring  # attempt to evaluate basestring
    def is_str(s):
        return isinstance(s, basestring)
except NameError:
    def is_str(s):
        return isinstance(s, str)

flatten = itertools.chain.from_iterable
import itertools
import gzip

def get_pymol_name(file_name):
    """Gets pymol name from file name"""
    name = os.path.basename(file_name)
    name = os.path.splitext(name)[0]
    #make it work with pdb.gz
    if name.endswith('.pdb'):
        name = name[:-4]
    name = cmd.get_legal_name(name)
    return name

def resolve_object_or_file_name(file_or_name):
    """Decides if the name is file path or an object name. 
    If it is a file path it loads the object and returns the name."""
    if os.path.exists(file_or_name):
        name = get_pymol_name(file_or_name)
        cmd.load(file_or_name, name)
    else:
        name = file_or_name
    return name


def read_json(filename):
    import json
    with opengz(filename, 'rt') as _file:
        return json.load(_file)


def get_ss(sel_str):
    """returns the secondary structure of selection as a string"""
    stored.ss_array = []
    cmd.iterate(f"{sel_str} and (name CA)", "stored.ss_array.append(ss)")
    return "".join(stored.ss_array)


def get_selection_property(sel_str, property="resi"):
    """returns the secondary structure of selection as a string"""
    stored.seli_array = []
    cmd.iterate(f"{sel_str} and (name CA)", f"stored.seli_array.append({property})")
    return stored.seli_array


get_sp = get_selection_property

def get_helices(sele="all"):
    helices = list()
    prevss = "nope"
    for atom in cmd.get_model(sele + " and name CA").atom:
        # print(atom.ss, atom.resi, atom.resi_number, atom.name)
        if atom.ss == "H":
            if atom.ss != prevss:
                helices.append(list())
            helices[-1].append(atom.resi)
        prevss = atom.ss
    return helices


def sel_first_helix(n, sele="all", onlyh=True, sele_name="sele"):
    helices = get_helices(sele)
    if int(onlyh):
        hresi = itertools.chain(*helices[: int(n)])
        sele_new = "resi " + "+".join(r for r in hresi)
    else:
        sele_new = "resi 0-%s" % helices[int(n) - 1][-1]

    sele_new = sele + " AND (" + sele_new + ")"
    cmd.select(sele_name, selection=sele_new)
    return sele_new


def sel_last_helix(n, sele="all", onlyh=True, sele_name="sele"):
    helices = get_helices(sele)
    if int(onlyh):
        hresi = itertools.chain(*helices[-int(n) :])
        sele_new = "resi " + "+".join(r for r in hresi)
    else:
        sele_new = "resi %s-99999" % helices[-int(n)][0]

    sele_new = sele + " AND (" + sele_new + ")"
    cmd.select(sele_name, selection=sele_new)
    return sele_new


def sel_terminal_helix(n, sele="all", onlyh=True, sele_name="sele"):
    if int(n) > 0:
        return sel_first_helix(int(n), sele=sele, onlyh=onlyh, sele_name=sele_name)
    else:
        return sel_last_helix(-int(n), sele=sele, onlyh=onlyh, sele_name=sele_name)



def sel_helix_from_to(from_n_1, to_n_1, sele="all", onlyh=True, sele_name=None):
    """Selects from_n to to_n helix, one based. Includes both from and to helices.
    Negative indices start counting from the end."""
    helices = get_helices(sele)
    N_hel = len(helices)
   

    if from_n_1>0:
        from_n=from_n_1-1
    else:
        from_n=N_hel+from_n_1

    if to_n_1>0:
        to_n=to_n_1-1
    else:
        to_n=N_hel+to_n_1
   
    if from_n>to_n:
        from_n,to_n=to_n,from_n
    
    #print(from_n,to_n)
    if int(onlyh):
        hresi = itertools.chain(*helices[from_n:to_n+1])
        sele_new = "resi " + "+".join(r for r in hresi)
    else:
        #from the first residue of the first helix, to the last residue of the last helix
        sele_new = f"resi {helices[from_n][0]}-{helices[to_n][-1]}"  

    sele_new = sele + " AND (" + sele_new + ")"
    if sele_name is not None: 
        cmd.select(sele_name, selection=sele_new)
    return sele_new

cmd.extend("sel_first_helix", sel_first_helix)
cmd.extend("sel_last_helix", sel_last_helix)
cmd.extend("sel_terminal_helix", sel_terminal_helix)
cmd.extend("sel_helix_from_to", sel_helix_from_to)

def clash_check_CA(selA, selB, distance=4.0, sele_name=None):
    """Returns the number of CA atoms in selA that lie closer than radii to any CA atom in selB"""
    ##
    selA = f"({selA}) and name CA"
    selB = f"({selB}) and name CA"

    sel_str = f"(({selA}) within {distance} of ({selB}))"
    model = cmd.get_model(sel_str)
    if not sele_name is None:
        cmd.select(sele_name, selection=sel_str)
    return model.nAtom

def get_distance2(c1, c2):
    return ((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)

def get_distance(c1, c2):
    return ((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2) ** 0.5


def get_alignment_map(from_sel, to_sel, max_distance=1.8):
    """Finds atoms in to_sel that are close to from_sel"""
    from_sel = f"({from_sel}) and name CA"
    to_sel = f"({to_sel}) and name CA"
    mapping = {}
    distances = {}
    from_model = cmd.get_model(from_sel)

    for at in from_model.atom:
        to_model = cmd.get_model(
            f"({to_sel}) within {max_distance} of ({from_sel} and chain {at.chain} and resi {at.resi})"
        )
        # print(f"({to_sel}) within {max_distance} of ({from_sel} and chain {at.chain} and resi {at.resi})")
        # print(to_model.nAtom)
        if to_model.nAtom > 1:
            #print(f"WARNING: more than one atom ({to_model.nAtom}) within {from_sel} and chain {at.chain} and resi {at.resi}. Choosing the closest")
            dist_pairs = []
            # sort the distances
            for to_atom in to_model.atom:
                dist_pairs.append( (get_distance2(to_atom.coord, at.coord), to_atom) )
            sorted_dist = sorted(dist_pairs, key = lambda x:x[0])
            to_model_dist, to_model_at = sorted_dist[0]
            to_model_dist = to_model_dist**0.5   
        elif to_model.nAtom == 1:
            to_model_at = to_model.atom[0]
            to_model_dist = get_distance(at.coord, to_model_at.coord)
        else:
            to_model_at = None
            to_model_dist = None
        
        ch_res_id = f"{at.chain}_{at.resi}"    
        if to_model_at: 
            mapping[ch_res_id] = f"{to_model_at.chain}_{to_model_at.resi}"
            distances[ch_res_id] = to_model_dist
        else:
            mapping[ch_res_id] = None
            distances[ch_res_id] = None
    return mapping, distances


def color_by_selector_array(name, selector, color="red"):
    if isinstance(selector, str):
        selector = selector.split(",")
    for n, i in enumerate(selector, start=1):
        if int(i) > 0:
            cmd.color(color, f"{name} and resi {n}")


def get_pymol_selector_from_rosetta_selector_str(sel_str):
    """Given a list of pdb_ids 10_A, or 10A print a pymol selector"""
    ids = sel_str.split(",")
    result = []
    # TODO group by chains?
    for id_ in ids:
        chain = id_[-1:]
        resi = id_[:-1]
        result.append(f"(resi {resi} and chain {chain})")
    result = "(" + " or ".join(result) + ")"
    return result
    
# get_pymol_selector_from_rosetta_selector_str('158B,160B,161B,162B,164B,165B,167B,168B,169B,171B,172B,194B,195B,196B,197B,198B,199B,200B,201B,202B,203B,204B,205B,206B,208B,209B')


def color_by_pdb_id_list(name, pdb_ids, color="red"):
    """Given a list of pdb_ids 10_A, or 10A print a pymol selector"""
    pass


cmd.extend("color_by_selar", color_by_selector_array)


def opengz(filename, mode):
    """Return either the system open or the gzip open """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def grep_file(file_name, expression):
    """Searches file and returns lines matching regular expression"""
    result = []
    
    with opengz(file_name, 'rt') as file:
        for line in file:
            #print(line)
            if re.search(expression, line): 
                result.append(line.strip())
    return result

def grep_lines(lines, expression):
    """Returns lines matchin expression"""
    result = []
    for line in lines:
        #print(line)
        if re.search(expression, line): 
            result.append(line.strip())
    return result


def read_file(filename):
    """Reads a file. Supports gzip files if ending is .gz"""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as myfile:
            data = myfile.read()
    else:
        with open(filename, 'rt') as myfile:
            data = myfile.read()
    return data  

#grep_file('out/repacked/ALAF05_8x/ALAF05_8x.pdb', r"^.*_pymol_selection")

def load_with_labels(pdb_path, print_cmd=False, prefix=''):
    """Loads a PDB and it's selections"""
    cmd_str=f"load {pdb_path}"
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
    import_labels_from_pdb(pdb_path, print_cmd=print_cmd, prefix=prefix)

load_wl = load_with_labels
cmd.extend("load_wl", load_wl)


def import_labels_from_pdb(pdb_path, print_cmd=False, ignore_res_seq_labels=True, load_labels=None, prefix='', object_name=None):
    if object_name is None:
        object_name = get_pymol_name(pdb_path)
    labels = grep_file(pdb_path, r"^REMARK PDBinfo-LABEL: ")
    #clean up first part
    labels = [lb.replace('REMARK PDBinfo-LABEL:','').strip() for lb in labels]
    # get all the different labels (split off the number and take the rest)
    if load_labels is None:
        unique_labels = sorted(set(flatten([lb.split()[1:] for lb in labels])))
    else:
        unique_labels = load_labels
    
    for label in unique_labels:
        if not label.startswith('res__') or not ignore_res_seq_labels:
            apply_label(labels, label, print_cmd=print_cmd, prefix=prefix, object_name=object_name)
   
apply_labels = import_labels_from_pdb
cmd.extend("apply_labels", apply_labels)
import_labels = import_labels_from_pdb
cmd.extend("import_labels", import_labels)

def apply_label(label_lines, label, print_cmd=False, prefix='', object_name=None):
    #must have spaces to avoid picking up similar labels
    choosen = grep_lines(label_lines, " "+label+"( |$)")
    res_num = [ch.split()[0] for ch in choosen]
    #print(label)
    #print(choosen)
    res_num = "+".join(res_num)
    if object_name is None:
        object_str = {}
    else:
        object_str = f"and {object_name}"
    cmd_str = f"select {prefix+label}, resi {res_num} {object_str}"
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)

def apply_read_selections(sel_str, print_cmd=False):
    """Applies read selection to pymol"""
    name, cmd_str = sel_str.split(' ', maxsplit=1)
    #The selections are named in newer versions of Rosetta
    #name = name.replace('_pymol_selection', "")
    #cmd_str = cmd_str.replace('rosetta_sele', name)
    print(name)
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
        
def import_selections_from_pdb(pdb_path, print_cmd=False):
    """Loads selection strings from the PDB file"""
    seles = grep_file(pdb_path, r"^.*_pymol_selection")
    for sel in seles:
        apply_read_selections(sel, print_cmd=print_cmd)

import_selections = import_selections_from_pdb
cmd.extend("import_selections", import_labels)

def load_with_selections(pdb_path, print_cmd=False):
    """Loads a PDB and it's selections"""
    cmd_str=f"load {pdb_path}"
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
    import_selections_from_pdb(pdb_path, print_cmd=print_cmd)

load_ws = load_with_selections
cmd.extend("load_ws", load_with_selections)

def apply_unsat_group(unsat_groups, object_name='all', print_cmd=False):
    """Parses an unsat group from BuriedUnsatHbonds (name: \n Unsatisfied HEAVY polar atom at residue 51: HIS  ND1 \n ... )"""
    lines = unsat_groups.split('\n')
    #name is the first line, skip the last colon
    name = lines.pop(0)[:-2]
    #print(name)
    selections = []
    for line in lines:
        line = line.strip()
        if line == '': continue #skip empty lines
        #split by spaces and take the last three fields
        spl = line.split()
        resnum = spl[-3][:-1] #skip the colum
        resname = spl[-2]
        atom   = spl[-1]
        
        sele = f"{resname}`{resnum}/{atom}"
        #print(sele)
        selections.append(sele)
    selections_str = " or ".join(selections)
    selections_str = f"(({selections_str}) and {object_name})"
    cmd_str = f"select {name}_unsats_{object_name}, {selections_str}"
    
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)


find_unsat_sections = re.compile(r"^BuriedUnsatHbonds\s(.*?)\s\sall_heavy_atom_unsats", re.MULTILINE | re.DOTALL )

def import_unsats_from_pdb(pdb_path, object_name=None, print_cmd=False):
    """Greps the PDB file for unsat blocks and turns them into pymol selections"""
    file_str = read_file(pdb_path)
    
    if object_name is None:
        object_name = get_pymol_name(pdb_path)
    unsat_groups = re.findall(find_unsat_sections, file_str)
    
    for unsat_group in unsat_groups:
        apply_unsat_group(unsat_group, object_name=object_name, print_cmd=print_cmd)

load_unsats = import_unsats_from_pdb
cmd.extend("load_unsats", load_unsats)
import_unsats = import_unsats_from_pdb
cmd.extend("import_unsats", import_unsats)

def import_per_residue_metric(pdb_path_or_lines, metric, object_name=None, print_cmd=False):
    """Loads a per_residue metric to the residues"""
    if is_str(pdb_path_or_lines):
        per_res_metrics = grep_file(pdb_path_or_lines, f"^{metric}_")
        if object_name is None:
            object_name = get_pymol_name(pdb_path_or_lines)
    else:
        per_res_metrics = grep_lines(pdb_path_or_lines, f"^{metric}_")
        if object_name is None:
            assert False, "Object name is not set when passing lines"

    for res in per_res_metrics:
        #strip the prefix, so we are left with a key-value pair
        res = res[len(metric)+1:]
        resi, val = res.split(' ', maxsplit=1)
        if print_cmd:
            print(f'alter {object_name} and resi {resi}, p.{metric}={val}')
        else:
            cmd.alter(f'{object_name} and resi {resi}', f"p.{metric}={val}")



def import_per_residue_metrics(pdb_path, metrics, object_name=None, print_cmd=False):
    """Loads a per_residue metric to the residues"""

    if object_name is None:
        object_name = get_pymol_name(pdb_path)
    
    if is_str(metrics):
        metrics = metrics.split()

    lines = read_file(pdb_path).split('\n')

    for metric in metrics:
        import_per_residue_metric(lines, metric, object_name)
    

import_per_res_metrics = import_per_residue_metrics
cmd.extend("import_per_res_metrics", import_per_res_metrics)
iprms = import_per_res_metrics 
cmd.extend("iprms", iprms)

def import_scores_from_pdb(pdb_path, object_name=None, print_cmd=False, 
        fields="total fa_rep fa_atr fa_sol fa_elec fa_dun_rot"):
    """
    Import the energy score table into pymol custom fields

    Parameters
    ----------
    pdb_path : str
        path to PDB
    object_name : str, optional
        the name of the object to load into, by default None
    print_cmd : bool, optional
        only print the command, by default False
    fields : str, optional
        [description], by default "total fa_rep fa_atr fa_sol fa_elec fa_dun_rot"
    """
    file_str = read_file(pdb_path)
    if object_name is None:
        object_name = get_pymol_name(pdb_path)

    #exrtract the energies table
    file_str=read_file(pdb_path)
    #Take everything between the tags. Must skip first line
    file_str=file_str.split('#BEGIN_POSE_ENERGIES_TABLE')[-1].split('#END_POSE_ENERGIES_TABLE')[0]
    #remove first line (#Begin...)
    file_str=file_str.split('\n', maxsplit=1)[-1]

    import csv, io
    reader = csv.DictReader(io.StringIO(file_str),  delimiter=' ', quotechar='"')

    if fields=="all" or fields is None:
        fields = reader.fieldnames
    else:
        #enable a string seperated with spaces
        if is_str(fields):
            fields=fields.split(' ')
        #check that the reqired fields are present
        fields_set = set(fields)
        reader_set = set(reader.fieldnames)
        assert fields_set.issubset(reader_set), f"Fields {fields_set.difference(reader_set)} are not in the PDB energy table. \
            \n fields    :{fields} \
            \n cvs_header:{reader.fieldnames}"

    fields_str = ""
    for field in fields:
        fields_str = fields_str + f"p.{field}=0;"
    if print_cmd:
        print(f'alter {object_name}, {fields_str}')
    else:
        cmd.alter(f'{object_name}', fields_str)

    for row in reader:
        label=row['label'] #label should be in the form GLU_18
        if label in ['weights', 'pose']: continue
        resi = label.split("_")[-1]
        fields_str = ""
        #in symetric structures a lot of rows will just be zero. Skip loading those
        all_zero = True
        for field in fields:
            fields_str = fields_str + f"p.{field}={row[field]};"
            all_zero = all_zero and (abs(float(row[field]))<1e-6)

        if not all_zero:
            if print_cmd:
                print(f'alter {object_name} and resi {resi}, {fields_str}')
            else:
                cmd.alter(f'{object_name} and resi {resi}', fields_str)

load_scores = import_scores_from_pdb
cmd.extend("load_scores", load_scores)
import_scores = import_scores_from_pdb
cmd.extend("import_scores", import_scores)

def color_by_score(field, range="-5 0 5", colors="green white red", selection="visible"):
    if is_str(range):
        range=range.split(' ')
        print(range)
    if is_str(colors):
       colors=colors.split(' ') 
    cmd.do(f'pseudoatom pOrig, pos=(0,0,0)')
    cmd.do(f'ramp_new proximityRamp, pOrig, selection=none, range={range}, color={colors}')
    cmd.do(f'spectrum p.{field}, {" ".join(colors)}, minimum={range[0]}, maximum={range[-1]}')

cmd.extend("color_by_score", color_by_score)
cmd.extend("cbs", color_by_score)


def load_error_prediction_data(npz_or_json_file):
    '''
        Metric types are:
    lddt, pe40, pe20, pe10, pe05
    '''
    import numpy as np
    import json
    if npz_or_json_file.endswith('.json') or npz_or_json_file.endswith('.json.gz') or \
       npz_or_json_file.endswith('.errpred') or npz_or_json_file.endswith('.errpred.gz'):
        data = read_json(npz_or_json_file) 
        
        #convert to np arrays
        for key in data.keys():
            data[key] = np.array(data[key])
        return data 
    
   
    dat = np.load(npz_or_json_file)
    
    lddt = dat["lddt"] 
    esto = dat["estogram"]
    res = {
        'lddt': lddt,
        'pe40': 1-np.mean(np.sum(esto[:,:,:4], axis=-1) + np.sum(esto[:,:,11:], axis=-1), axis=-1),
        'pe20': 1-np.mean(np.sum(esto[:,:,:5], axis=-1) + np.sum(esto[:,:,10:], axis=-1), axis=-1),
        'pe10': 1-np.mean(np.sum(esto[:,:,:6], axis=-1) + np.sum(esto[:,:,9:], axis=-1), axis=-1),
        'pe05': 1-np.mean(np.sum(esto[:,:,:7], axis=-1) + np.sum(esto[:,:,8:], axis=-1), axis=-1),
    }

    return res
    
def apply_metric_from_npz(npz_file, metric_type='lddt', sel_str='all', print_cmd=False):
    '''Imports an array of values from the npz file
    Metric types are:
    lddt, pe40, pe20, pe10, pe05
    '''
    import numpy as np
    dat = load_error_prediction_data(npz_file)
    res = dat[metric_type]

    model = cmd.get_model(f"{sel_str} and name CA")
    
    assert (len(res)==model.nAtom)

    for n, atom in enumerate(model.atom):
        cmd_str = f"alter {sel_str} and resi {atom.resi}, b='{res[n]}'"
        if not print_cmd:
            cmd.do(cmd_str)
        else:
            print(cmd_str)
apply_erp = apply_metric_from_npz
cmd.extend("apply_erp", apply_metric_from_npz)
#import is the recommended naming
cmd.extend("import_erp", apply_metric_from_npz)
import_metric_from_npz = apply_metric_from_npz
import_erp = apply_metric_from_npz

def load_with_error_metrics(pdb_path, metric_type='lddt', selection=None, npz_path=None, print_cmd=False):
    """Loads a PDB and the error metrics"""
    
    if npz_path is None:
        npz_path = pdb_path.replace('.pdb', '.npz')

    if selection is None:
        selection = get_pymol_name(pdb_path)     

    cmd_str=f"load {pdb_path}"
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
    apply_metric_from_npz(npz_path, metric_type=metric_type, sel_str=selection, print_cmd=print_cmd)
        
load_erp = load_with_error_metrics
cmd.extend("load_erp", load_with_error_metrics)

r"""
dela
run Z:\projects\truncator\truncator\pymol_utils.py
load_erp N1_C3_DR64_HS57_cryoat_mALb8_cutT1_BA__c_term__o7__r34_0001_0004_0001_asym.pdb, metric_type=lddt
spectrum b, rainbow, minimum=0.2, maximum=1
print(get_alignment_map("/ZCON_1__numH3__from-22.38__to07.23__grAB-CD//A","DHR08_trim"))
"""


def load_rosetta_pdb(pdb, scores=True, labels=True, unsats=True, per_res_metrics=[], color_by="score", object=None):
    """
    Loads a rosetta specific pdb with all the additional info.

    Parameters
    ----------
    pdb : str
        pdb file 
    scores : bool, optional
        load scores, by default True
    labels : bool, optional
        load labels as selections, by default True
    unsats : bool, optional
        load unsats as dots, by default True
    """
    if object is None:
        object = get_pymol_name(pdb)

    cmd.load(pdb, object=object)
    if labels:
        import_labels(pdb, object_name=object)
    if unsats:
        import_unsats(pdb, object_name=object)
        cmd.show('dots', f'*unsats_{object}*')
        #cmd.show('dots', f'*unsats_{object}*')
    if scores:
        import_scores(pdb, object_name=object)

    if per_res_metrics:
        import_per_residue_metrics(pdb, metrics=per_res_metrics, object_name=object)

    if color_by:
        if color_by=='score':
            color_by=='total'
        color_by_score(color_by)
    

cmd.extend("load_rosetta_pdb", load_rosetta_pdb)
lrp = load_rosetta_pdb
cmd.extend("lrp", lrp)


# taken from http://www.protein.osaka-u.ac.jp/rcsfp/supracryst/suzuki/jpxtal/Katsutani/InterfaceResidues.py
from pymol import stored
 
def interface_residues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    """
    interfaceResidues -- finds 'interface' residues between two chains in a complex.
 
    PARAMS
        cmpx
            The complex containing cA and cB
 
        cA
            The first chain in which we search for residues at an interface
            with cB
 
        cB
            The second chain in which we search for residues at an interface
            with cA
 
        cutoff
            The difference in area OVER which residues are considered
            interface residues.  Residues whose dASA from the complex to
            a single chain is greater than this cutoff are kept.  Zero
            keeps all residues.
 
        selName
            The name of the selection to return.
 
    RETURNS
        * A selection of interface residues is created and named
            depending on what you passed into selName
        * An array of values is returned where each value is:
            ( modelName, residueNumber, dASA )
 
    NOTES
        If you have two chains that are not from the same PDB that you want
        to complex together, use the create command like:
            create myComplex, pdb1WithChainA or pdb2withChainX
        then pass myComplex to this script like:
            interfaceResidues myComlpex, c. A, c. X
 
        This script calculates the area of the complex as a whole.  Then,
        it separates the two chains that you pass in through the arguments
        cA and cB, alone.  Once it has this, it calculates the difference
        and any residues ABOVE the cutoff are called interface residues.
 
    AUTHOR:
        Jason Vertrees, 2009.		
    """
    # Save user's settings, before setting dot_solvent
    oldDS = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)
 
    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName+"1"
    chA, chB = "chA", "chB"
 
    # operate on a new object & turn off the original
    cmd.create(tempC, cmpx)
    cmd.disable(cmpx)
 
    # remove cruft and inrrelevant chains
    cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
 
    # get the area of the complete complex
    cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    cmd.alter(tempC, 'q=b')
 
    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    cmd.extract(chA, tempC + " and (" + cA + ")")
    cmd.extract(chB, tempC + " and (" + cB + ")")
    cmd.get_area(chA, load_b=1)
    cmd.get_area(chB, load_b=1)
 
    # update the chain-only objects w/the difference
    cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
 
    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')
 
    cmd.enable(cmpx)
    cmd.select(selName1, None)
    for (model,resi,diff) in stored.r:
        key=resi+"-"+model
        if abs(diff)>=float(cutoff):
            if key in seen: continue
            else: seen.append(key)
            rVal.append( (model,resi,diff) )
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))
 
    # this is how you transfer a selection to another object.
    cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    cmd.delete(selName1)
    cmd.delete(chA)
    cmd.delete(chB)
    cmd.delete(tempC)
    # show the selection
    cmd.enable(selName)
 
    # reset users settings
    cmd.set("dot_solvent", oldDS)
 
    return rVal
 
cmd.extend("interface_residues", interface_residues)


def get_chainbreaks(objname_or_file, cutoff_A=2, cmd=None, delete_all_before_load=False):
    """Returns the resi where the chain break is. The next residue is in a different place 
    N --chain_break-- N+1  .Returns N."""
    import numpy as np
    if cmd is None:
        import pymol
        cmd = pymol.cmd

    if os.path.exists(objname_or_file):
        objname = get_pymol_name(objname_or_file)
        if delete_all_before_load:
            cmd.delete('all')
        cmd.load(objname_or_file, object=objname)
    else:
        objname = objname_or_file 

    C_atoms = cmd.get_coords(objname+" and name C", 1)
    N_atoms = cmd.get_coords(objname+" and name N", 1)
    #subastract the C from the N of the next residue
    distances = np.sum((C_atoms[:-1]-N_atoms[1:])**2, axis=1)
    #len(distances), min(distances), np.mean(distances) ,max(distances)
    breaks = distances > cutoff_A**2
    return breaks.nonzero()[0]+1


aa_3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

def get_seq_alignment_map(from_sele, to_sele, aln_name=None, cmd=pymol.cmd, cycles_align=5, cutoff_align_RMSD=2.5):
    """Gets the sequnce mapping from calling align. Returns a tuple of CA atoms"""
    model1_name = cmd.get_object_list(from_sele)
    assert len(model1_name)==1
    model1_name = model1_name[0]
    model2_name = cmd.get_object_list(to_sele)
    assert len(model2_name)==1   
    model2_name = model2_name[0]

    if aln_name is None:
        aln_name = f"aln__{model1_name}__{model2_name}"
    
    cmd.align(f'({to_sele}) and name CA', f'({from_sele}) and name CA', cutoff=cutoff_align_RMSD, cycles=cycles_align, object=aln_name)

    #cmd.super(f'({to_sele}) and name CA', f'({from_sele}) and name CA', object=aln_name)
    raw_aln  = cmd.get_raw_alignment(aln_name)
    res = [] 
    for idx1, idx2 in raw_aln:
        from_at = cmd.get_model(f"index {idx1[1]} and {idx1[0]}").atom[0]
        from_at.parent_model = idx1[0]
        to_at = cmd.get_model(f"index {idx2[1]} and {idx2[0]}").atom[0]
        to_at.parent_model = idx2[0]

        if   model1_name == idx1[0] and  model2_name == idx2[0]:
            res.append((from_at, to_at))
        elif model2_name == idx1[0] and  model1_name == idx2[0]:
            res.append((to_at, from_at))
        else:   
            raise "Models names don't match"
    


    return res

def print_seq_alignement_map(smap):
    """Prints a tuple of atoms as alignment map"""
    model1 = smap[0][0].parent_model 
    model2 = smap[0][1].parent_model
    
    print(f'{model1} -> {model2}')
    for at1, at2 in smap:
        print(f'{at1.chain} {at1.resi_number:3d} ({aa_3to1[at1.resn]}) -> {at2.chain} {at2.resi_number:3d} ({aa_3to1[at2.resn]})')

def get_seq_aln_map_differences(smap):
    """Returns only the tuples that have a diffrent resisue name (resn) """
    res=[]
    for at1, at2 in smap:
        if at1.resn != at2.resn:
            res.append((at1, at2))
    return res


def get_selections_from_seq_aln_map(smap):
    """Returns the from and to selections strings (based on atom index and name)"""

    #each selection can only be from one model, so only look at first atom
    #TODO rewrite using get_object_list
    model1 = smap[0][0].parent_model 
    model2 = smap[0][1].parent_model

    sel1_idx = []
    sel2_idx = []
    for at1, at2 in smap:
        sel1_idx.append(f'({model1} and chain {at1.chain} and resi {at1.resi})')
        sel2_idx.append(f'({model2} and chain {at2.chain} and resi {at2.resi})')
    sel1_idx_str = ' or '.join(sel1_idx)
    sel2_idx_str = ' or '.join(sel2_idx)
    

    #TODO: REWRITE SO IT'S GROUPED BY CHAIN
    sel1 = f'byresi( {sel1_idx_str} )' 
    sel2 = f'byresi( {sel2_idx_str} )' 
    
    return sel1, sel2
    
def filter_seq_aln_map_on_selection(smap, sel, smap_index=0):
    """Only keep those alignments that are both in sel1 and sel2"""
    sel_model = cmd.get_model(f"{sel} and name CA")
    sel_indexes = set([atom.index for atom in sel_model.atom])
    res = []
    for link in smap:    
        if link[smap_index].index in sel_indexes:
            res.append(link)
    
    return res

def show_align_diff(sel1, sel2, col1='same', col2='same', make_sel=True, print_sel=True, cmd=pymol.cmd):
    """Finds the diffrences in sequence (after alignment) and colors them if color is given"""
    smap = get_seq_alignment_map(sel1, sel2) 
    #print_seq_alignement_map(smap)

    smap = get_seq_aln_map_differences(smap)
    

    if len(smap) == 0:
        print('NO DIFFERENCES')
        return False
    diff_sel1, diff_sel2 = get_selections_from_seq_aln_map(smap)
    
    model1 = smap[0][0].parent_model 
    model2 = smap[0][1].parent_model
    
    if make_sel:
        cmd.select(f'{model1}_diff_{model2}', diff_sel1)
        cmd.select(f'{model2}_diff_{model1}', diff_sel2)
        
    if not col1 is None:
        if col1 != 'same':
            cmd.color(col1, diff_sel1)
        cmd.show('licorice', diff_sel1)
        cmd.hide('licorice', f" ({diff_sel1}) and (h. and (e. c extend 1))")
    if not col2 is None:
        if col2 != 'same':
            cmd.color(col2, diff_sel2)
        cmd.show('licorice', diff_sel2)
        cmd.hide('licorice', f" ({diff_sel2}) and (h. and (e. c extend 1))")
    
    cmd.do('util.cnc')

    if print_sel:
        print_seq_alignement_map(smap)        
    return True

cmd.extend("show_align_diff", show_align_diff)

def mutate_residue(sele, target_3resname, cmd=pymol.cmd):
    """Mutates the sele to target_3resname"""
    #Open wizard
    cmd.do("wizard mutagenesis; refresh_wizard") 
    
    cmd.select('sele', f'byresi ({sele})')
    cmd.get_wizard().do_select('''sele''')
    cmd.get_wizard().do_state(8)
    cmd.get_wizard().set_mode(target_3resname)
    cmd.get_wizard().do_state(1)
    
    
    cmd.get_wizard().apply()
    #TODO: Copy chi angles here
    #close wizard
    cmd.set_wizard()



cmd.extend("mutate_residue", mutate_residue)


def transfer_residue(from_sele, to_sele, copy_rotamer=False):
    """Mutates residue 'to_sele' to the same type as 'from_sele'
    from_sele and to_sele must only be one residues.
    
    TODO: Copy rotamer angles, better error message
    """
    from_res = cmd.get_model(from_sele +' and name CA')
    assert  len(from_res.atom)==1, 'only one from res allowed but got more'
    from_res = from_res.atom[0]
    
    to_res = cmd.get_model(to_sele +' and name CA')
    assert  len(to_res.atom)==1, 'only one from res allowed but got more'
    to_res = to_res.atom[0]
    
    #print(from_res.resn)
    mutate_residue(to_sele, from_res.resn)

cmd.extend("transfer_residue", transfer_residue)

def pymol_display(cmd=pymol.cmd, ray=False):
    """Displays the pymol session in a Jupyter notebook"""
    
    from IPython.display import Image
    from tempfile import mktemp
    import os
    image_name = mktemp()+'.png'
    cmd.png(image_name, ray=ray)
    im = Image(filename=image_name)
    os.remove(image_name)
    return im


def get_rmsd(ob_file1, ob_file2, align=False, cmd=None, sub_sel='chain A and name CA', cleanup=True)->int:
    """
    gets rmsd without aligning

    Parameters
    ----------
    ob_file1 : file or obj name
        file or obj name
    ob_file2 : file or obj name
        file or obj name
    align : bool, optional
        Align before calculating the RMSD, by default False
    cmd :
        pymol cmd module, by default None
    sub_sel : str, optional
        substring to select, by default 'chain A and name CA'
    cleanup : bool, optional
        delete models after comparing, by default True

    Returns
    -------
    int
        RMSD
    """    
    """Returns RMSD."""
    if cmd is None:
        import pymol
        cmd = pymol.cmd
    
    if os.path.exists(ob_file1):
        ob1 = get_pymol_name(ob_file1)
        cmd.load(ob_file1, object=ob1)
    else:
        ob1 = ob_file1 
        
    if os.path.exists(ob_file2):
        ob2 = get_pymol_name(ob_file2)
        cmd.load(ob_file2, object=ob2)
    else:
        ob2 = ob_file2    
        
    #return None
    #cmd.remove('not chain A')
    if align:
        rmsd =  cmd.rms(f'{ob1} and {sub_sel}', f'{ob2} and {sub_sel}', cutoff=100, cycles=0)
    else:
        #cmd_str = f"rms_cur {ob1} and {sub_sel}, {ob2} and {sub_sel}, cutoff=100, cycles=10"
        #cmd.do(cmd_str)
        rmsd = cmd.rms_cur(f'{ob1} and {sub_sel}', f'{ob2} and {sub_sel}', cutoff=100, cycles=0)
       

    if cleanup:
        cmd.delete(ob1)
        cmd.delete(ob2)
        
    return rmsd


def get_hfuse_rmsd(hfuse_pdb, block1, block2, align_block1=True, align_block2=True, cleanup=True, cmd=None, hfuse_extra_sel=''):
    """returns RMSD of hfuse. Returns a dict of values. The H-fuse must have the overlap label set"""
    if cmd is None:
        import pymol
        cmd = pymol.cmd
    import numpy as np
    
    hfuse = resolve_object_or_file_name(hfuse_pdb)
    #overlap selection gets made
    import_labels_from_pdb(hfuse_pdb, load_labels=['overlap']) 

    block1 = resolve_object_or_file_name(block1)
    block2 = resolve_object_or_file_name(block2)

    #print(block1, block2)
    
    if align_block1:
        cmd.align(block1, hfuse+' '+hfuse_extra_sel, gap=-20.0, extend=-0.2, max_gap=500)
    if align_block2:
        cmd.align(block2, hfuse+' '+hfuse_extra_sel, gap=-20.0, extend=-0.2, max_gap=500)
    
    #cmd.super(block1, hfuse)
    #cmd.super(block2, hfuse)

    alm1, dist1 = get_alignment_map('overlap', block1, max_distance=200)
    distn1 = np.array(list(dist1.values()))
    rmsd1 = np.sqrt(np.average(distn1**2))

    alm2, dist2 = get_alignment_map('overlap', block2, max_distance=200)
    distn2 = np.array(list(dist2.values()))
    rmsd2 = np.sqrt(np.average(distn2**2))

    #print(dist1, dist2)
    #take the bigger rmsd.
    if rmsd1>rmsd2:
        rmsd=rmsd1
        distn=distn1
    else:
        rmsd=rmsd2
        distn=distn2

    result = {
        'overlap_rmsd' :rmsd,
        'overlap_n_all':len(distn),
        'overlap_n_gt_05':sum(distn>0.5),
        'overlap_n_gt_10':sum(distn>1.0),
        'overlap_n_gt_15':sum(distn>1.5),
        'overlap_n_gt_20':sum(distn>2.0),
    }

    if cleanup:
        cmd.delete(hfuse)
        cmd.delete(block1)
        cmd.delete(block2)
        
    return result

def fuse_on_helix_overlap(obj1, obj2, fusion_term, obj1_chain='A', obj2_chain='A', 
                            obj1_overlap_hel=(1,2), obj2_overlap_hel=(1,2), cmd=None, verbose=False):
    """
    Used to fuse to structures with a helical overlap. Especially useful for designed helical repeat proteins (DHRs).
    Obj1 is modified in place.
    Obj2 is modifed and deleted.

    Returns:
    --------
    Modified Obj1 and a dictionary with overlap metrics.
       'overlap_rmsd': 0.2, #overlap of the selections, 
       'aligned_atom_num': 274, 
       'aligned_res_num': 55, 
       'NC_length_A': 1.6 # Length of the new NC bond in Angstrem

    Parameters
    ----------
    obj1 : str or file name
        Existing object name or file to load.
        This is the target that will be fused on.
    obj2 : str or file name
        This gets deleted.
    fusion_term : str
        Which terminal of Obj1 to fuse to (n_term, c_term)
    obj1_chain : str, optional
        Chain to fuse on, by default 'A'
    obj2_chain : str, optional
        Chain that is fused, by default 'B'
    obj1_overlap_hel : tuple, optional (1,2).
        Which helices to use as overlap.
        If c_term, these are counted form the end.
        One based counting.
    obj2_overlap_hel : tuple, optional
        If c_term, these are counted form the end.
        One based counting.
    """
    if cmd is None:
        import pymol
        cmd = pymol.cmd

    obj1 = resolve_object_or_file_name(obj1)
    obj2 = resolve_object_or_file_name(obj2)

    res = {}

    obj1_sel = f'{obj1} and chain {obj1_chain}'
    obj2_sel = f'{obj2} and chain {obj2_chain}'

    if fusion_term=='c_term':
        obj1_overlap_hel=-obj1_overlap_hel[0],-obj1_overlap_hel[1]
    else:
        obj2_overlap_hel=-obj2_overlap_hel[0],-obj2_overlap_hel[1]

    obj1_overlap = f'{sel_helix_from_to(obj1_overlap_hel[0], obj1_overlap_hel[1], obj1_sel, onlyh=False)} AND backbone '
    obj2_overlap = f'{sel_helix_from_to(obj2_overlap_hel[0], obj2_overlap_hel[1], obj2_sel, onlyh=False)} AND backbone'
    #print(obj1_overlap)
    #print(obj2_overlap)
    cmd.color('tv_red', obj1_overlap)
    cmd.color('red', obj2_overlap)
    if verbose:
        print('obj1_overlap:', obj1_overlap)
        print('obj2_overlap:', obj2_overlap)
    
    aln_res = cmd.align(obj2_overlap, obj1_overlap, cycles=0)

    res['overlap_rmsd'] = aln_res[0]
    res['aligned_atom_num'] = aln_res[1]
    res['aligned_res_num'] = aln_res[-1]
    if verbose:
        print(f'ALN_RES: {aln_res}')

    ### remove duplicate chains
    cmd.alter('obj2', 'segi="fused"')

    obj1_chains = cmd.get_chains(obj1)
    obj2_chains = cmd.get_chains(obj2)

    #print(obj1_chains, obj2_chains)
    # remove duplicate chains (except the one that is being fused, that will get changed anyway)
    # the chains get a lower case for now
    # TODO make the chain switching algo more robust 
    for an_obj2_chain in obj2_chains:
        if (an_obj2_chain in obj1_chains) and (an_obj2_chain != obj2_chain):
            new_chain_id = an_obj2_chain.lower()
            if verbose:
                print(f"Changing {an_obj2_chain} to {new_chain_id}")
                print(f'{obj2} and chain {an_obj2_chain}', f'chain="{new_chain_id}"')
            cmd.alter(f'{obj2} and chain {an_obj2_chain}', f'chain="{new_chain_id}"')

    #delete redundant residues
    #TODO Get the smallest distance between atoms and do the crossover there.
    if fusion_term=='c_term':
        obj1_overlap_mdl = cmd.get_model(f'({obj1_overlap}) AND name C')
        last_obj1 = list(obj1_overlap_mdl.atom)[-1]
        cmd.remove(f"{obj1} and chain {last_obj1.chain} and segi '{last_obj1.segi}' and resi {last_obj1.resi_number+1}-" )

        obj2_overlap_mdl = cmd.get_model(f'({obj2_overlap}) AND name N')
        last_obj2 = list(obj2_overlap_mdl.atom)[-1]
        cmd.remove(f"{obj2} and chain {last_obj2.chain} and segi '{last_obj2.segi}' and resi 1-{last_obj2.resi_number}" )

        #change the fusion numbering
        cmd.alter(f'{obj2} and chain {last_obj2.chain}',  f'resi=int(resi)-{last_obj2.resi_number}+{last_obj1.resi_number}')
        
        #change the fusion chain abd segi
        cmd.alter(f'{obj2} and chain {last_obj2.chain}',  f'segi="{last_obj1.segi}"')
        cmd.alter(f'{obj2} and chain {last_obj2.chain}',  f'chain="{last_obj1.chain}"')

        C_atom = f"{obj1} and chain {last_obj1.chain} and segi '{last_obj1.segi}' and resi {last_obj1.resi_number} and name C"
                                    #not a mistake             #not a mistake                   #not a mistake    
        N_atom = f"{obj2} and chain {last_obj1.chain} and segi '{last_obj1.segi}' and resi {last_obj1.resi_number+1} and name N"
        if verbose:
            print("N_ATOM: ", N_atom)
            print("C_ATOM: ", C_atom)
        cmd.show('spheres', N_atom)
        cmd.show('spheres', C_atom)
        
        res['NC_length_A'] = cmd.get_distance(N_atom, C_atom)
        cmd.fuse(N_atom, C_atom, move=0, mode=0)
        #cmd.bond(N_atom, C_atom)

    else:
        raise "n_term fusion not yet implemented. So sorry." 

    cmd.delete('obj2')
    res['obj1'] = obj1
    res['obj2'] = obj2
    return res