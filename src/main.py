# $ sage main.py 

if __name__ == "__main__":
    from func_test import Test
    time1_ps,time2_ps=Test("Proposed","Square",1)
    time1_po,time2_po=Test("Proposed","One",1)
    time1_es,time2_es=Test("Existing","Square",1)
    time1_eo,time2_eo=Test("Existing","One",1)
    time1_ps,time2_ps=Test("Proposed","Square",2)
    time1_po,time2_po=Test("Proposed","One",2)
    time1_es,time2_es=Test("Existing","Square",2)
    time1_eo,time2_eo=Test("Existing","One",2)
