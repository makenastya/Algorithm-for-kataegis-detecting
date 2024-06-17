def get_kataegis(path):
    """Принимает на вход путь до файла .maf, возвращает репорт о найденных катаегисах в виде DataFrame."""
    maf = pd.read_csv(path, sep='\t')
    if 'type' in maf.columns:
        df = maf[['type', 'Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']]
    else:
        df = maf[['Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']]
    df.loc[df['Variant_Type']. isin(['INS', 'DEL', 'SNP'])]
    g = df.groupby(['Tumor_Sample_Barcode'])
    if 'type' in df.columns:
        rep = pd.DataFrame(columns = ['type', 'SampleName', 'chr', 'start', 'end', 'width', 'totalVariants'])
    else:
        rep = pd.DataFrame(columns = ['SampleName', 'chr', 'start', 'end', 'width', 'totalVariants'])
    
    for i in g:
        sample_data = i[1]
        name = i[0][0]
        chrom = sample_data.groupby(['Chromosome'])
        for chr_data in chrom:
            chr_data_sort = chr_data[1].sort_values(by='Start_Position')
            if len(chr_data_sort) >= 6:
                m = 1
                while m < len(chr_data_sort):
                    k = 1
                    if 'type' in rep.columns:
                        type = chr_data_sort.iloc[m - 1]['type']
                    start = chr_data_sort.iloc[m - 1]['Start_Position']
                    while m < len(chr_data_sort) and chr_data_sort.iloc[m]['Start_Position'] - chr_data_sort.iloc[m - 1]['End_Position'] <= 1000:
                        k += 1
                        m += 1
                    end = chr_data_sort.iloc[m - 1]['End_Position']
                    if k >= 6:
                        if 'type' in rep.columns:
                            rep.loc[len(rep.index)] = [type, name, chr_data[0][0], start, end, end - start + 1, k]
                        else:
                            rep.loc[len(rep.index)] = [name, chr_data[0][0], start, end, end - start + 1, k]
                    m += 1
    return(rep)
