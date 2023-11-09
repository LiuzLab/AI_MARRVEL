select DISTINCT phm.acc_num, phm.gene_sym, phm.phen_id, con.code
        from hgmd_phenbase.hgmd_mutation phm
        join hgmd_phenbase.phenotype_concept phc on phm.phen_id = phc.phen_id 
        join hgmd_phenbase.concept con on con.cui = phc.cui  where con.sab = "HPO" limit 10
        INTO OUTFILE '/var/lib/mysql-files/HGMD_phen.tsv'
        FIELDS TERMINATED BY '\t';

