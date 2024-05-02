def get_census_codes(minority):
    
    assert minority in {'Asian', 'Black', 'Hispanic', 'Native'}
    
    # Census codes available at http://starr.tamu.edu/files/2013/01/Census-Codes.pdf
    if minority == 'Asian':
        codes = ['P0030006'] # Asian alone
        codes += ['P0030013','P0030017','P0030020','P0030023','P0030024'] # Asian (among 2 races)
        codes += ['P0030028','P0030031','P0030034','P0030035','P0030037','P0030040','P0030041','P0030043','P0030044','P0030046'] # 3
        codes += ['P0030048','P0030051','P0030052','P0030054','P0030055','P0030057','P0030058','P0030059','P0030061','P0030062'] # 4
        codes += ['P0030064','P0030065','P0030067','P0030068','P0030069'] # 5
        codes += ['P0030071'] # 6
    elif minority == 'Black':
        codes = ['P0030004'] # Black or African American alone 
        codes += ['P0030011','P0030016','P0030017','P0030018','P0030019'] # Black or African American (among 2 races)
        codes += ['P0030027','P0030028','P0030029','P0030030','P0030037','P0030038','P0030039','P0030040','P0030041','P0030042'] # 3
        codes += ['P0030048','P0030049','P0030050','P0030051','P0030052','P0030053','P0030058','P0030059','P0030060','P0030061'] # 4
        codes += ['P0030064','P0030065', 'P0030066','P0030067','P0030069'] # 5
        codes += ['P0030071'] # 6
    elif minority == 'Hispanic':
        codes = ['P0040002'] # Hispanic or Latino VAP
    elif minority == 'Native':
        codes = ['P0030005'] # American Indian and Alaska Native alone
        codes += ['P0030012','P0030016','P0030020','P0030021','P0030022'] # American Indian and Alaska Native (among 2 races)
        codes += ['P0030027','P0030031','P0030032','P0030033','P0030037','P0030038','P0030039','P0030043','P0030044','P0030045'] # 3
        codes += ['P0030048','P0030049','P0030050','P0030054','P0030055','P0030056','P0030058','P0030059','P0030060','P0030062'] # 4
        codes += ['P0030064','P0030065','P0030066','P0030068','P0030069'] # 5
        codes += ['P0030071'] # 6
            
    return codes
