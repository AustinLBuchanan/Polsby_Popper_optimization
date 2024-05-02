from carve import carving_heuristic
from complete import recursive_bipartition_heuristic
from cleanup import mip_local_search
from metrics import report_metrics

# All three phases in one function!
# 1. Carve
# 2. Complete
# 3. Cleanup

def gingles(G, minority, mvap_lower=0.5, mvap_upper=1.0, mvap_excess_penalty=0.001, 
            polsby_popper_cutoff=0.125, complement_compactness=False, deviation_penalty=0.0, h=1, preserve_splits=True):
    
    # Phase 1: carve Gingles districts
    gingles_districts = carving_heuristic(G, minority, mvap_lower=mvap_lower, mvap_upper=mvap_upper, 
                                          mvap_excess_penalty=mvap_excess_penalty, polsby_popper_cutoff=polsby_popper_cutoff, 
                                          complement_compactness=complement_compactness, deviation_penalty=deviation_penalty)
    print(f"Finished carving off {len(gingles_districts)} Gingles districts.")

    # Phase 2: complete the plan
    complete_plan = recursive_bipartition_heuristic(G, gingles_districts=gingles_districts)
    print("\nFor completed plan:")
    report_metrics(G, complete_plan, minority)
    
    # Phase 3: cleanup the plan 
    cleaned_plan = complete_plan.copy()
    for hops in range(1,h+1):
        cleaned_plan = mip_local_search(G, cleaned_plan, minority=minority, h=hops, preserve_splits=preserve_splits)
        print(f"\nFor cleaned plan ({hops} hops):")
        report_metrics(G, cleaned_plan, minority)
    
    return cleaned_plan
