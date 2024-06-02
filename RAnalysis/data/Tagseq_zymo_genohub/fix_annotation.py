input_file = "Pver_genome_assembly_v1.0.gff3"
output_file = "Pver_genome_assembly_v1.0_modified.gff3"
log_file = "annotation_repair.log"

# Notes:
# .gff3 file is separated by gene blocks - i.e. until you reach a new gene, it is safe to assume that all current
#   annotations are related to the last gene you saw.
# I don't want to assume that there can be only 1 mrna per gene. There is also a similar assumption as to the gene
#   assumption; until you reach a new mRNA/gene annotation, all annotations are related to the last mRNA you saw.
# The problem we are trying to solve is that for some reason these non-mRNA and non-gene annotations have the wrong
#   parent ID. They are, in some cases, listing their own ID as the parentID.

# Method to the madness:
# - Scan the file for a 'mRNA' tag
# - Until you reach another mRNA tag:
#   - read in every annotation
#   - check that its parent ID is set to the mRNA ID you last saw
#   - If it is: 
#       spit into the next file without comment.
#   - If it isn't:
#       spit into the next file with the correct parent ID, take note.

# Unless otherwise stated, every line is copied into the output file exactly as-is.


# Helper function.
# Takes in the string of the entirety of the final column, splits it into the individual attributes.
# Then extracts the ID and Parent attributes, and returns them as a tuple.
def get_id_pid(datcol) :
    # sub-columns within the final column are separated by a semicolon
    cols = datcol.split(';')

    # Get only the attributes 'ID' and 'Parent'
    i_id = None
    i_parent = None
    for col in cols :
        # Every attribute is stored as 'name=value', so by splitting on '=' we can get back the name of the attribute
        #   and the value of the attribute.
        (attr, val) = col.split("=")
        if attr == "ID" :
            i_id = val
        elif attr == "Parent" :
            i_parent = val

    # If we could not find values for this feature's ID and its parent's ID, throw an exception. Otherwise, return
    #   the values we found.
    if i_id == None or i_parent == None :
        raise ValueError("Could not find both ID and parent ID in data: " + datcol)
    else :
        return(i_id, i_parent)

# Helper function.
# Takens in the string of the entirety of the annotation column, and what to change the parent ID to.
# Creates an array that stores all of the attr=value strings of all attributes in the column, except for the parent
#   attribute. Once it reaches the parent attribute, it replaces it with an attr=value string containing the new_parent
#   id for the value.
# Then returns a string that contains all of these attr=value strings combined with semicolons separating them.
def adjust_pid(datcol, new_parent) :
    cols = datcol.split(';')
    output = []
    for col in cols :
        (attr, val) = col.split("=")
        if(attr == "Parent") :
            output.append("Parent=" + new_parent)
        else :
            output.append(col)
    return(";".join(output) + "\n")

# Helper function.
# Takes in an array of data; the first index is the correctly formatted data column that follows the gff standard.
#   Every following index contains a space-separated pair of attribute and its value that needs to be added to the 
#   correctly-formatted column.
# Example:
# INPUT     ["ID=cds.Pver_g25526.t1;Parent=Pver_g25526.t1", "5_prime_partial true", "3_prime_partial true"]
# OUTPUT    "ID=cds.Pver_g25526.t1;Parent=Pver_g25526.t1;5_prime_partial=true;3_prime_partial=true"
def fix_extra_columns(datcols) :
    # Get the correctly formatted first index
    output = datcols[0]
    # For every following index, split by space and add to the output string
    for col in datcols[1:] :
        split = col.split(" ")
        output = output + ";" + split[0] + "=" + split[1]
    return output

cur_mrna_id = "NA"
cur_id = None
cur_parent = None

with open(input_file, 'r') as inp :
    with open(output_file, 'w') as out :
        with open(log_file, 'w') as log:
            for lineno, line in enumerate(inp, 1):
                # If it starts with a '#', ignore it and spit it into the output file. These are just comments.
                if line[0] == "#" :
                    out.write(line)
                    continue

                columns = line.strip().split("\t")

                dat_type = columns[2]
                
                # If there are more columns than expected, that means we have the additional 5_prime_partial or
                #   3_prime_partial tags incorrectly appended. Fix the attribute column to store this data, and then
                #   get rid of the extra columns. Log that the change has been made.
                if(len(columns) > 9) :
                    columns[8] = fix_extra_columns(columns[8:])
                    columns = columns[:9]
                    line = "\t".join(columns) + "\n"
                    log.write("Modified line: " + str(lineno) + "\tFixed misplaced attributes.\n")

                attributes = columns[8]


                # if the 3rd column is not gene or mRNA -- and we have not yet seen an mrna -- throw an error.
                #   This is only because of the assumptions I made up in the notes section; I want to know if they're not
                #   valid! If I can't assume the block-like structure where every exon, CDS, etc is related to the mRNA
                #   annotation proceeding it, we flat out can't repair the file because we have no idea what the real
                #   mRNA parent is.
                if cur_mrna_id == "NA" and dat_type != "gene" and dat_type != "mRNA":
                    raise AssertionError("Assumption of file structure is incorrect. Repair cannot proceed.")

                # if the type is gene, just paste it into the next file without anything fancy.
                if dat_type == "gene" :
                    out.write(line)
                # If the type name is mRNA, update our stored mRNA name and paste it into the output file without any
                #   changes.
                elif dat_type == "mRNA":
                    cur_mrna_id, cur_parent = get_id_pid(attributes)
                    out.write(line)
                # If it is of any other type, check to make sure that the parent ID is set correctly.
                # If not, log it in our log file & write to the output file an adjusted line.
                else :
                    cur_id, cur_parent = get_id_pid(attributes)
                    if(cur_parent != cur_mrna_id) : #the parent is supposed to be the id of the mrna we saw last
                        log.write("Modified line: " + str(lineno) + "\tTo: " + cur_mrna_id + "\tFrom: " + cur_parent + "\n")
                        tmp = columns[:8]
                        tmp.append(adjust_pid(columns[8], cur_mrna_id))

                        out.write("\t".join(tmp))
                    else :
                        out.write(line)