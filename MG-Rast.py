"""
Author: Mike Aaldering
Company: NIOO Knaw
Email: M.Aaldering2@nioo.knaw.nl
Private Email: mike.aaldering098@gmail.com

This program downloads species and metadata from MG-Rast database and stores it to a .CSV file.
"""
import collections, csv, ctypes, ijson, json, re , six, time, urllib2
import pandas as pd
from dateutil import parser
from ete2 import NCBITaxa
from incf.countryutils import transformations


def main():

    only_standardize = raw_input('Only standardize type: "n", \ndownload data type: "d" \n')
    CSVfile_name = "MG-Rast"                           # CSV file where your metadata and species are getting saved
    only_species_filename = CSVfile_name + " Only Species.txt"
    only_metadata_filename = CSVfile_name + " Only Metadata.txt"
    fileinterestingspecies = "InterestingSpecies.txt"   # FILE where species of interest are saved
    material = "soil,wetland,sediment,peat"             # Material of interest
    limit = "50"                                        # Limit of samples downloading per time
    relative = True                                     # Make species relative
    standardize = True                                  # If you want to standardize the metadata set it to True
    fieldname_order_filename = "Data_Order.txt"
    fieldname_order = order_fielnames(fieldname_order_filename)
    interestingspecieslist = get_species_of_interest(fileinterestingspecies, fileinterestingspecies)

    if only_standardize == "n":
        print "Standardizing data"
        join_data(only_species_filename, only_metadata_filename, standardize, fieldname_order, CSVfile_name + " Only standardized")
    elif only_standardize == "d":
        url = make_url(limit, material)
        total_count = get_total_count(url, limit)
        all_metadata = []
        all_species = []
        howfar = 0
        while (howfar <= (int(total_count)-int(limit))):
            print "Downloading: ", howfar, "-", int(howfar) + int(limit), " of total: ", total_count
            response = urllib2.urlopen(url)
            data = json.load(response)
            for sample in data["data"]:
                ############GETTING METADATA#########
                metadata_dict = flatten_json(sample)
                metadata_dict = no_new_line(metadata_dict)
                all_metadata.append(metadata_dict)

                ############GETTING SPECIES#########
                species_dict = get_species(metadata_dict["id"], relative, interestingspecieslist)
                all_species.append(species_dict)


                ############WRITING SPECIES AND METADATA TO SEPARATE FILE#################
                write_csv(all_species, only_species_filename, "0")
                write_csv(all_metadata, only_metadata_filename, "NA")

                ######standardize DATA AND JOIN DATA#######
                join_data(only_species_filename, only_metadata_filename, standardize, fieldname_order, CSVfile_name)

                print howfar, " / ", total_count, metadata_dict["id"]
                howfar = howfar + 1


            url = data["next"]
    else:
        print only_standardize, "is not validate, try it again"
        main()
    print "done!"

def join_data(only_species_filename, only_metadata_filename, standardize, fieldname_order, CSVfile_name):
    ##########JOIN METADATA WITH SPECIES###################
    df_species = pd.read_csv(only_species_filename, error_bad_lines=False, sep='\t')
    df_metadata = pd.read_csv(only_metadata_filename, error_bad_lines=False, sep='\t')
    MetaSpecies = pd.concat([df_metadata, df_species], axis=1, join_axes=[df_metadata.index])
    MetaSpecies_list_with_dict = MetaSpecies.T.to_dict().values()
    all_data_list_with_dict = make_list_witch_dict_without_nan(MetaSpecies_list_with_dict)

    ######standardize DATA#######
    if standardize == True:
        final_all_data_list_with_dict = []
        for row in all_data_list_with_dict:
            row_dict = standardize_data(row)
            final_all_data_list_with_dict.append(row_dict)
    else:
        final_all_data_list_with_dict = all_data_list_with_dict

    ##########WRITE FILE WITH ALL DATA###############
    CleanDataFrame = pd.DataFrame(final_all_data_list_with_dict)
    header_list = list(CleanDataFrame.columns.values)
    data_order = order_data(header_list, fieldname_order)
    CleanDataFrame = CleanDataFrame[data_order]
    CleanDataFrame = CleanDataFrame.fillna(value="NA")

    all_data_header = get_level(list(CleanDataFrame.columns.values))
    DataFrame2 = pd.DataFrame(all_data_header)
    DataFrame2 = DataFrame2.set_index(["Orig"])
    DataFrame2 = DataFrame2.T
    DataFrame2 = DataFrame2.append(CleanDataFrame)
    DataFrame2.to_csv(CSVfile_name + ".txt", sep='\t', encoding="utf-8")


def make_list_witch_dict_without_nan(list_with_dict):
    all_data = []
    for row in list_with_dict:
        row_dict = {}
        for data in row:
            if str(row[data]) != "nan":
                row_dict[data] = row[data]
                #print boem, iets[boem]
        all_data.append(row_dict)
    #print all_data
    return all_data

def write_csv(data, filename, fill):
    DataFrame = pd.DataFrame(data)
    DataFrame = DataFrame.fillna(value=fill)
    DataFrame.to_csv(filename, sep='\t', encoding="utf-8", index=False)


def get_level(header):
    all_data_header = []
    for species in header:
        species_one = []
        species_one.append(species)
        ncbi = NCBITaxa()

        taxids = ncbi.get_name_translator(species_one)
        lineage_list = []
        if bool(taxids) == True:
            lineage_list = ncbi.get_lineage(taxids[species][0])
        name_list = ncbi.translate_to_names(lineage_list)
        header_dict = {}
        header_dict["Orig"] = species
        for idx, item in enumerate(name_list):
            header_dict[idx] = item
        all_data_header.append(header_dict)
    # print all_data_header
    return all_data_header


def order_data(header_list, fieldname_order):
    new_list = []
    for item in fieldname_order:
        if item in header_list:
            new_list.append(item)
            # fieldname_order = fieldname_order.remove(item)
    for deze in header_list:
        if deze not in new_list:
            new_list.append(deze)
    return new_list


def order_fielnames(filename):
    fieldnames = []
    with open(filename, "r") as f:
        for row in f:
            dataname = row.rstrip('\n')
            # print dataname
            if dataname not in fieldnames:
                fieldnames.append(dataname)
    return fieldnames


def no_new_line(metadata_dict):
    for datapoint in metadata_dict:
        standardize = True
        try:
            # Deleting new lines in data
            if "\n" in str(metadata_dict[datapoint]):
                metadata_dict[datapoint.lstrip()] = re.sub("\n", '', str(metadata_dict[datapoint]))
            if "\r" in str(metadata_dict[datapoint]):
                metadata_dict[datapoint.lstrip()] = re.sub("\r", '', str(metadata_dict[datapoint]))
            # Renaming all the column headers
            if "metadata sample data" in datapoint:
                new_key = re.sub("metadata sample data", '', datapoint)
                new_key = new_key + "_sample"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata library data" in datapoint:
                new_key = re.sub("metadata library data", '', datapoint)
                new_key = new_key + "_library"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata env_package data" in datapoint:
                new_key = re.sub("metadata env_package data", '', datapoint)
                new_key = new_key + "_env_package"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata project data" in datapoint:
                new_key = re.sub("metadata project data", '', datapoint)
                new_key = new_key + "_project"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata env_package data" in datapoint:
                new_key = re.sub("metadata env_package data", '', datapoint)
                new_key = new_key + "_env_package"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata env_package" in datapoint:
                new_key = re.sub("metadata env_package", '', datapoint)
                new_key = new_key + "_env_package"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata library data" in datapoint:
                new_key = re.sub("metadata library data", '', datapoint)
                new_key = new_key + "_library"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata library" in datapoint:
                new_key = re.sub("metadata library", '', datapoint)
                new_key = new_key + "_library"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata project data" in datapoint:
                new_key = re.sub("metadata project data", '', datapoint)
                new_key = new_key + "_project"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata project" in datapoint:
                new_key = re.sub("metadata project", '', datapoint)
                new_key = new_key + "_project"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata sample data" in datapoint:
                new_key = re.sub("metadata sample data", '', datapoint)
                new_key = new_key + "_sample"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata sample data" in datapoint:
                new_key = re.sub("metadata sample data", '', datapoint)
                new_key = new_key + "_sample"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
            elif "metadata sample" in datapoint:
                new_key = re.sub("metadata sample", '', datapoint)
                new_key = new_key + "_sample"
                metadata_dict[new_key.lstrip()] = metadata_dict.pop(datapoint)
        except UnicodeEncodeError:
            pass

            # print "We can't encode: ", datapoint, metadata_dict[datapoint], " this is not a big problem."

    return metadata_dict



def get_metadata_of_interest(file):
    with open(file, "r") as f:
        metadatalist = []
        for metadata in f:
            nonewlinedata = metadata.rstrip('\n')
            metadatalist.append(nonewlinedata)
        # print specieslist
        return metadatalist


def get_species_of_interest(species, file):
    with open(file, "r") as f:
        specieslist = []
        for species in f:
            nonewlinespecies = species.rstrip('\n')
            specieslist.append(nonewlinespecies)
            firstword = nonewlinespecies.split(' ')
            specieslist.append(firstword[0])
        specieslist = list(set(specieslist))
        # print specieslist
        return specieslist



def standardize_data(metadata_dict):
    #country
    if "country_sample" in metadata_dict:
        if "French" in str(metadata_dict["country_sample"]):
            metadata_dict["country_sample"] = "France"
        if "United Kingdom" in str(metadata_dict["country_sample"]):
            metadata_dict["country_sample"] = "UK"
        if "United Status of America" in str(metadata_dict["country_sample"]):
            metadata_dict["country_sample"] = "USA"
        if "Viet Nam" in str(metadata_dict["country_sample"]):
            metadata_dict["country_sample"] = "Vietnam"
        if " " in str(metadata_dict["country_sample"]):
            metadata_dict["country_sample"] = str(metadata_dict["country_sample"]).replace(" ", "_")

    #env_package_sample
    if "env_package_sample" in metadata_dict:
        if "biofilm" in str(metadata_dict["env_package_sample"]):
            metadata_dict["env_package_sample"] = "biofilm"
    #pH
    if "ph_sample" in metadata_dict:
        metadata_dict["pH_std"] = metadata_dict["ph_sample"]
    if "ph_env_package" in metadata_dict:
        metadata_dict["pH_std"] = metadata_dict["ph_env_package"]

    if "pH_std" in metadata_dict:
        if "; pH" in str(metadata_dict["pH_std"]):
            value = metadata_dict["pH_std"]
            value = value.encode('utf-8')
            value = value.replace("; pH", "")
            metadata_dict["pH_std"] = value
        if "-" in str(metadata_dict["pH_std"]):
            ph = metadata_dict["pH_std"]
            ph = ph.split("-")
            ph = ((float(ph[1])-float(ph[0]))/2) + float(ph[0])
            metadata_dict["pH_std"] = ph
    # Nitrogen
    if "tot_n_env_package" in metadata_dict:
        if "%" in str(metadata_dict["tot_n_env_package"]):
            metadata_dict["N_percentage_std"] = metadata_dict["tot_n_env_package"].replace("%", "")
        if "mg/kg nitrate" in str(metadata_dict["tot_n_env_package"]):
            value = metadata_dict["tot_n_env_package"]
            value = value.replace("mg/kg nitrate", "")
            metadata_dict["N_percentage_std"] = (float(value) * 14) / 620000


    # Carbon
    if "_tot_org_c_env_package" in metadata_dict:
        value = metadata_dict["_tot_org_c_env_package"]
        if "%" in str(value):
            metadata_dict["C_percentage_std"] = value.replace("%", "")
        else:
            metadata_dict["C_percentage_std"] = value
    if "org_carb_env_package" in metadata_dict:
        value = metadata_dict["org_carb_env_package"]
        if "(%wt)" in str(value):
            metadata_dict["C_percentage_std"] = value.replace("(%wt)", "")
        else:
            metadata_dict["C_percentage_std"] = value
    if "tot_org_carb_env_package" in metadata_dict:
        if metadata_dict["id_project"] == "mgp6200":
            metadata_dict["C_percentage_std"] = metadata_dict["tot_org_carb_env_package"]
    # Description
    if "metadata project data project_description" in metadata_dict:
        if metadata_dict["id"] == "mgm4441091.3" or metadata_dict["id"] == "mgm4441619.3" or metadata_dict[
            "id"] == "mgm4441620.3" or metadata_dict["id"] == "mgm4441656.4":
            metadata_dict["metadata project data project_description"] = "NA"
    # Temperature
    if "temperature_sample" in metadata_dict:
        if "; Celsius" in str(metadata_dict["temperature_sample"]):
            value = metadata_dict["temperature_sample"]
            value = value.replace("; Celsius", "")
            metadata_dict["Temp_std"] = value
        elif "-" in str(metadata_dict["temperature_sample"]):
            value = str(metadata_dict["temperature_sample"])
            if value[0] == "-":
                metadata_dict["Temp_std"] = metadata_dict["temperature_sample"]
            else:
                # value = value.encode('utf-8')
                value = value.split("-")
                value = ((float(value[1]) - float(value[0])) / 2) + float(value[0])
                metadata_dict["Temp_std"] = value
        else:
            metadata_dict["Temp_std"] = metadata_dict["temperature_sample"]
    # sequence methode
    if "seq_make_library" in metadata_dict:
        if "454" in str(metadata_dict["seq_make_library"]):
            metadata_dict["seq_std"] = "pyrosequencing"
        if "Roche" in str(metadata_dict["seq_make_library"]):
            metadata_dict["seq_std"] = "pyrosequencing"
        if "HiSeq" in str(metadata_dict["seq_make_library"]):
            metadata_dict["seq_std"] = "Illumina"
        else:
            metadata_dict["seq_std"] = metadata_dict["seq_make_library"]
    # alitude
    if "metadata sample data altitude" in metadata_dict:
        value = metadata_dict["metadata sample data altitude"]
        key = "altitude_std"
        metadata_dict[key] = value
    # country
    if "metadata sample data country" in metadata_dict:
        value = str(metadata_dict["metadata sample data country"])
        key = "metadata sample data country"
        if value == "USA":
            value = "United States of America"
            metadata_dict[key] = value
            key = "continent_std"
            value = transformations.cn_to_ctn(value)
            metadata_dict[key] = value
        elif value == "Pacific Ocean":
            value = "Pacific_Ocean"
            metadata_dict[key] = value
            metadata_dict["continent_std"] = value
        elif value == "England":
            value = "Europe"
            key = "continent_std"
            metadata_dict[key] = value
        elif value == "United Kingdom":
            metadata_dict[key] = "England"
            value = "Europe"
            key = "continent_std"
            metadata_dict[key] = value
        elif value == "Viet Nam":
            metadata_dict[key] = "Viet Nam"
            value = "Asia"
            key = "continent_std"
            metadata_dict[key] = value
        elif value == "UK":
            metadata_dict[key] = "England"
            value = "Europe"
            key = "continent_std"
            metadata_dict[key] = value
        elif value == "South Korea":
            metadata_dict[key] = "South Korea"
            value = "Asia"
            key = "continent_std"
            metadata_dict[key] = value
        elif value == "French Republic":
            value = "France"
            metadata_dict[key] = value
            key = "continent_std"
            metadata_dict[key] = ["Europe"]
        else:
            value = metadata_dict["metadata sample data country"]
            value = transformations.cn_to_ctn(value)
            key = "continent_std"
            metadata_dict[key] = value
    if "metadata sample data collection_date" in metadata_dict:
        value = metadata_dict["metadata sample data collection_date"]
        value.splitlines()
        try:
            d = parser.parse(value)
            value = d.strftime("%Y")
            key = "collection_year_std"
            metadata_dict[key] = value
        except ValueError:
            print "error"
    if "metadata library data seq_meth" in metadata_dict:
        value = str(metadata_dict["metadata library data seq_meth"])
        if value == "454":
            # print value
            value = "pyrosequencing"
        key = "seq_meth_std"
        metadata_dict[key] = value
    return metadata_dict


def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def get_species(sampleid, relative, interestingspecieslist):
    url = "http://api.metagenomics.anl.gov/1/metagenome/" + sampleid + "?verbosity=stats"
    response = urllib2.urlopen(url)
    objects = ijson.items(response, '')
    for item in objects:
        species_list_dict = item["statistics"]["taxonomy"]["species"]
        species_dict = {}
        for item in species_list_dict:
            species_dict[item[0]] = item[1]

        if relative == True:
            total_reads = 0
            total_reads_interested = 0
            relative_species_dict = {}
            for species in species_dict:
                total_reads = total_reads + species_dict[species]
            interesting_species_dict = only_interesting_species(species_dict, interestingspecieslist)
            for species in interesting_species_dict:
                total_reads_interested = total_reads_interested + species_dict[species]

            for species in interesting_species_dict:
                relative = (float(interesting_species_dict[species]) / float(total_reads_interested))
                relative_species_dict[species] = relative
            som = 0
            for raar in relative_species_dict:
                som = som + float(relative_species_dict[raar])
            relative_species_dict["total_reads"] = total_reads
            relative_species_dict["total_reads_interested"] = total_reads_interested
            return relative_species_dict
        else:
            interesting_species_dict = only_interesting_species(species_dict, interestingspecieslist)
            return interesting_species_dict


def only_interesting_species(species_dict, interestingspecies_list):
    keys = []
    for key in species_dict:
        for species in interestingspecies_list:
            species_low = key.lower()
            interesting_low = species.lower()
            if interesting_low in species_low:
                keys.append(key)
    species_dict = {x: species_dict[x] for x in keys}

    to_remove = []
    for name in species_dict:
        if "Methylobacterium" in name:      #Methylobacterium is no methanotroph
            to_remove.append(name)
    for remove in to_remove:
        del species_dict[remove]

    return species_dict


def make_url(limit, material):
    url = "http://api.metagenomics.anl.gov/metagenome?sequence_type=WGS&metadata=" + material + "&limit=" + limit + "&verbosity=metadata"
    return url


def get_total_count(url, limit):
    splitter = "limit=" + str(limit)
    url = url.split(splitter)
    url = url[0] + "limit=1" + url[1]
    response = urllib2.urlopen(url)
    data = json.load(response)
    print "total samples downloading: ", data['total_count']
    return data['total_count']


def flatten_json(y):
    out = {}
    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + ' ')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + ' ')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


main()
