from number_of_districts import congressional_districts_2020
from gerrychain import Graph
from datetime import date
import networkx as nx
import gurobipy as gp
import os, sys, math, csv
import export, ordering, hess
import mip, mip_contiguity, mip_objective, mip_fixing, mip_local_search
datapath = '..\\districting-data-2020\\'

def main():
    
    # Check args
    args = sys.argv[1:]
    if args_okay(args):
        print("Using arguments:",args)
        (state, level, objective, contiguity) = args
    else:
        return
    
    # get today's date
    today = date.today()
    today_string = today.strftime("%Y_%b_%d") # Year_Month_Day, like 2019_Sept_16
    
    # export filepath
    export_filepath = '../results_' + today_string + '/'
    
    # export results to a row of the summary csv file
    summary_csv = export_filepath + 'summary.csv'
    my_fieldnames = ['state','level','objective','contiguity'] # arguments
    my_fieldnames += ['k','L','U','n','m'] # params
    my_fieldnames += ['B_size', 'B_time', 'hess_time', 'ls_time'] # max B and heuristic info
    my_fieldnames += ['MIP_obj','MIP_bound','MIP_time', 'MIP_status', 'MIP_nodes', 'callbacks', 'lazy_cuts'] # MIP info
    
    # if results directory and csv file doesn't exist yet, then create them
    if not os.path.exists(export_filepath):
        os.makedirs(export_filepath)
        with open(summary_csv,'w',newline='') as csvfile:   
            writer = csv.DictWriter(csvfile, fieldnames = my_fieldnames)
            writer.writeheader()
    
    # store experimental results in 'result' dictionary
    result = dict()
    (result['state'],result['level'],result['objective'],result['contiguity']) = args
    
    # Read graph G and identify number of districts k.
    filename = state + '_' + level + '.json'
    G = Graph.from_json( datapath + filename )
    if not nx.is_connected(G):
        print("ERROR: graph is disconnected. Aborting...")
        return
    DG = nx.DiGraph(G) # directed version of G
    DG._state = state
    DG._level = level
    for node in DG.nodes:
        DG.nodes[node]['TOTPOP'] = G.nodes[node]['P0010001']
    
    # What is ideal district population?
    DG._k = congressional_districts_2020[state]
    ideal = sum( DG.nodes[i]['TOTPOP'] for i in DG.nodes ) / DG._k 
    
    # Calculate lower and upper population bounds L and U
    #   under a 1% population deviation (+/-0.5%).
    dev = 0.01
    DG._L = math.ceil( ( 1 - dev / 2 ) * ideal )
    DG._U = math.floor( ( 1 + dev / 2 ) * ideal )
    print("Using L =",DG._L,"and U =",DG._U,"and k =",DG._k)
    
    (result['k'],result['L'],result['U']) = ( DG._k, DG._L, DG._U ) 
    (result['n'],result['m']) = ( DG.number_of_nodes(), G.number_of_edges() )
    
    # Build base MIP model
    m = mip.build_base_mip(DG)
    m._numCallbacks = 0
    m._numLazyCuts = 0
    m._callback = None
    
    # Add objective (and any related constraints)
    if objective == 'cut':
        mip_objective.add_cut_edges_objective(m, DG)
    elif objective == 'perim':
        mip_objective.add_perimeter_objective(m, DG)
    elif objective == 'invpp':
        mip_objective.add_inverse_Polsby_Popper_objective(m, DG)
    elif objective == 'avepp':
        mip_objective.add_average_Polsby_Popper_objective(m, DG)
    else:
        print("ERROR: this should not happen. This objective not supported:",objective)
    
    # Find max B set and vertex ordering,
    #   for use w/ symmetry handling and variable fixing
    (B, B_time) = ordering.solve_maxB_problem(DG)
    (result['B_size'], result['B_time']) = ( len(B), '{0:.2f}'.format(B_time) )
    DG._ordering = ordering.find_ordering(DG, B)
    
    # Add contiguity constraints
    if contiguity == 'lcut':
        m._DG = DG
        m.Params.LazyConstraints = 1
        m._callback = mip_contiguity.lcut_callback
    elif contiguity == 'scf':
        mip_contiguity.add_scf_constraints(m, DG)
    elif contiguity == 'shir':
        mip_contiguity.add_shir_constraints(m, DG)
    else:
        print("ERROR: this should not happen. These contiguity constraints not supported:",contiguity)
    
    # Solve Hess model (perhaps heuristically) to get a feasible solution
    if level == 'county':
        (hess_labeling, hess_time) = hess.solve_hess_model(DG)
    elif level == 'tract':
        (hess_labeling, hess_time) = hess.hess_heuristic(DG)
    else:
        hess_time = 0.00
        print("ERROR: warm start heuristic assumes level in {county,tract}")
    result['hess_time'] = '{0:.2f}'.format(hess_time)
    
    # MIP-based local search
    if hess_labeling:
        # convert hess_labeling into one that meets the partitioning orbitope restrictions
        orbitope_friendly_labeling = mip.get_orbitope_friendly_labeling(DG, hess_labeling)
        
        # Improve solution quality with MIP-based local search,
        #   but first add one-root-per-district constraints!
        #   They are implied by partitioning orbitope EF, but EF hasn't been added yet.
        root_constrs = m.addConstrs( gp.quicksum( m._r[i,j] for i in DG.nodes ) == 1 for j in range(DG._k) )
        (ls_labeling, ls_time) = mip_local_search.local_search(m, DG, orbitope_friendly_labeling, radius=1)
        m.remove(root_constrs)
        
        # Inject local search warm start
        mip.inject_warm_start(m, DG, ls_labeling)
        
        result['ls_time'] = '{0:.2f}'.format(ls_time)
    else:
        result['ls_time'] = 'n/a'
    
    # Symmetry handling (partitioning orbitope)
    mip.add_partitioning_orbitope_constraints(m, DG)
    
    # Add variable fixings
    mip_fixing.do_variable_fixing(m, DG)
    
    # Solve
    m.Params.TimeLimit = 600
    m.Params.MIPGap = 0.00
    m.Params.IntFeasTol = 1.e-9
    m.Params.FeasibilityTol = 1.e-9
    m.optimize(m._callback)
    
    # Solution reporting
    result['MIP_time'] = '{0:.2f}'.format(m.runtime)
    result['MIP_status'] = int(m.status)
    result['MIP_nodes'] = int(m.NodeCount)
    result['MIP_bound'] = m.objBound
    result['callbacks'] = m._numCallbacks
    result['lazy_cuts'] = m._numLazyCuts
    
    # report best solution found, if any
    if m.SolCount <= 0:
        result['MIP_obj'] = 'no_solution_found'
    else:
        result['MIP_obj'] = m.objVal
        labeling = { i : j for i in DG.nodes for j in range(DG._k) if m._x[i,j].x > 0.5 }
        
        # export districting map to png (image)
        export_filename = export_filepath + state + '_' + level + '_' + objective + '_' + contiguity
        export.export_to_png(DG, labeling, export_filename + '.png')

        # export districting plan to block assignment file (csv)
        export.export_to_baf(DG, labeling, export_filename + '.baf')
    
    export.append_dict_as_row( summary_csv, result, my_fieldnames )
    
    return
    
    
def args_okay(args):
    
    if len(args) != 4:
        print("Expecting 4 arguments: <state> <level> <objective> <contiguity>")
        return False

    state = args[0]
    level = args[1]
    objective = args[2]
    contiguity = args[3]
    
    state_args = [ key for key in congressional_districts_2020.keys() ]
    level_args = [ 'county', 'tract' ]
    objective_args = [ 'cut', 'perim', 'invpp', 'avepp' ]
    contiguity_args = [ 'scf', 'lcut', 'shir' ] 
    
    if state not in state_args:
        print(state,"is not a permitted state. Pick from",state_args)
        return False
    elif level not in level_args:
        print(level,"is not a permitted level. Pick from",level_args)
        return False
    elif objective not in objective_args:
        print(objective,"is not a permitted objective. Pick from",objective_args)
        return False
    elif contiguity not in contiguity_args:
        print(contiguity,"is not a permitted contiguity constraint. Pick from",contiguity_args)
        return False
    else:
        return True
   
        
if __name__ == "__main__":
    main()
    