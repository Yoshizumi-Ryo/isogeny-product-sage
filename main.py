# $ sage main.py 

if __name__ == "__main__":
    from func_test import Test
    time1_ps,time2_ps=Test("Proposed","Square")
    time1_po,time2_po=Test("Proposed","One")
    time1_es,time2_es=Test("Existing","Square")
    time1_eo,time2_eo=Test("Existing","One")
