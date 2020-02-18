using DataFrames, DataFramesMeta, TableView
#using CSV (no longer needed)

# This program processes Riskman sequence output and divides each initiator's
# sequences into groups according to the primary failure that causes core damage.

# test_sort4.jl does not output CSV files. Instead it is for data analysis
# purposes to help with the writeups of differences between RiskMan and FinPSA
# results.

# represents a basic event found in the Riskman sequence output
struct BasicEvent
    name::String
    description::String
    probability::AbstractFloat # occurrence probability (or frequency)
    isfailure
end

# I didn't end up using this struct to represent sequences.
# Instead I just made them rows in a dataframe.
struct Seq
    number::Int
    initiator::String
    frequency::AbstractFloat
    keyfailure::String # the name of the basic event
    events::Array{BasicEvent}
end

# function takes a line from a CSV file and splits it into an array of
# strings, one for each cell. In cells that contain commas and therefore have
# quotes in the CSV file, it removes the quotes
function unCSV(line::String)
    if ! ('"' in line)
        return split(line,',')
    end
    line = replace(line, "\"\"\"" => "\"")
    pieces = split(line,'"')
    if length(pieces) == 3
        return vcat(split(pieces[1],',')[1:end-1],
            pieces[2],
            split(pieces[3],',')[2:end])
    else
        println("Error in unCSV, line has multiple quoted cells.")
    end
end

function readfiles()
    lines = String[]
    # open up all six CSV files with the Riskman sequence results
    # and concatenate them into a big list of lines
    for i in 1:6
        counter = 0
        print("File ", i)
        filename = "./all_seqs$i.csv"
        open(filename) do file
            for line in eachline(file)
                push!(lines,line)
                counter += 1
            end
            print(" ", counter, " lines, ", length(lines), " total\n")
        end
        # frame = CSV.read(filename)
        #@show first(frame,3)
    end
    return lines
end

# This function goes through and divides the lines from the CSV files into 
# Riskman sequences and the split fractions (basic events) that constitute them.
# It also makes a list of all the initiating events used in these sequences.
function process_lines(lines)
    BEs = Dict() # unique basic events
    inits = Dict() # initiating events
    seqs = DataFrame(Seqnum=[], Init=[], Freq=[], Failure=[], Evts=[])
    seqnum = 0
    counter = 0
    failindex = -1
    evt_list = []
    for line in lines
        counter += 1
        if counter % 1000 == 0
            println(counter, " ", line)
        end
        pieces = unCSV(line)
        # throw out lines that don't contain an event
        #println(pieces[1], " ", pieces[2], " ", pieces[3])
        if !(tryparse(Float64, pieces[1]) isa Number)
            continue
        end
        be_name = pieces[3]
        if ! haskey(BEs, be_name)
            # a basic event we haven't seen before. Add it to the list.
            BEs[be_name] = BasicEvent(pieces[3], pieces[6]*pieces[7], tryparse(Float64, pieces[4]), nothing)
        end
        newnum = tryparse(Int64, pieces[1])
        if newnum != seqnum
            # make a new sequence (wait why did i have this sequence object?)
            #newseq = Seq(newnum, BEs[be_name], tryparse(Float64,pieces[10]), [], -1)
            seqnum = newnum
            evt_list = [] # I keep a pointer to evt_list after putting it in the df
                # then I can add events to the list later on as I find them.
            failindex = -1
            # add new sequence to the sequence dataframe
            push!(seqs, (seqnum, be_name, tryparse(Float64,pieces[10]), "", evt_list))
            # add the basic event to the initiating events list
            if ! haskey(inits, be_name)
                inits[be_name] = BEs[be_name]
            end
        else
            # add event to current sequence
            push!(evt_list, BEs[be_name])
            # also look for the key failure (the one with the smallest probability)
            if failindex < 1 || evt_list[end].probability < evt_list[failindex].probability
                failindex = length(evt_list)
                # Set the Failure column for the current sequence to the name of the BE
                # with the lowest failure probability
                seqs[(seqs[!,:Seqnum] .== seqnum),:Failure] .= evt_list[failindex].name
            end
        end
    end
    return BEs, inits, seqs
end

function group_by_init(seqs, inits)
    seqs_by_init = Dict()
    for init in keys(inits)
        # also sort by frequency, high to low
        seqs_by_init[init] = sort!(seqs[(seqs[!,:Init] .== init),:], [:Freq], rev=true)
    end
    return seqs_by_init
end

function group_by_failure(seqs_by_init)
    groupings = Dict()
    groups = Dict()
    println("Number of sequences for each initiator: ")
    for (k,v) in seqs_by_init
        println(k, " ", size(v,1)) # print number of sequences for initiator k
        groupings[k] = sort!(by(seqs_by_init[k], :Failure, :Freq => sum), :Freq_sum, rev=true)
        #showtable(seqs_by_init[k])
        #showtable(groupings[k])
        initfails = groupings[k][!,:Failure]
        groups[k] = []
        for f in initfails
            push!(groups[k],@where(seqs_by_init[k], :Failure .== f))
        end
    end
    return groups
end

# find if seq is a non-minimal version of any previous sequence in the group
# returns a string saying which sequence it's a non-minimal version of, and
# a list of the failure(s) that make it non-minimal (the extra failures)
function nonminimal(seq, group, groupnum)
    nonmin = ""
    verbose = false
    nonmin_events = []
    for seq2 in eachrow(group)
        if seq2 == seq
            break # can't be a non-minimal version of a sequence w/ lower frequency
        end
        all_repeated = true # start out assuming seq contains every event in seq2
        for evt in seq2[:Evts]
            if !(evt in seq[:Evts]) # if an event in seq2 turns out not to be in seq
                all_repeated = false # then seq is not a non-minimal version of seq2
                verbose && (nonmin = string("failed on ", evt))
                break
            end
        end
        if all_repeated
            nonmin = string(", Non-,Minimal (Seq ", seq2[:Seqnum], ")\n")
            nonmin_events = setdiff(seq[:Evts], seq2[:Evts])
            if seq[:Init] == "SMLOCA" && groupnum == 1
                verbose && println("Checking seq ", seq[:Seqnum], " against ", seq2[:Seqnum], " : ", nonmin)
            end
            break
        end
        if seq[:Init] == "SMLOCA" && groupnum == 1
            verbose && println("Checking seq ", seq[:Seqnum], " against ", seq2[:Seqnum], " : ", nonmin)
        end
    end
    if startswith(nonmin, "failed")
        nonmin = ""
    end
    return (nonmin, nonmin_events)
end

function find_repeats(group)
    all_evts = []
    for seq in eachrow(group)
        push!(all_evts, [be.name for be in seq[:Evts]])
    end
    # for a in all_evts
    #     println(a)
    # end
    repeats = intersect(all_evts...)
    return repeats
end

# This function is for outputting all the sorted sequences into CSV files,
# which can be imported into Excel. The stuff in the first line is all
# formatting hints for creating conditional formatting rules in Excel, so
# that for instance any split fraction that makes the sequence non-minimal 
# can be colored blue.
function print_to_csv(groups, BEs)
    # nonmins = ["CLATEC", "CLATED", "DPACC1", "OAMIS1", "CFCFB", "DCH", "TIVR1"]
    # Each initiator will have multiple sequence groups
    for (init,grouplist) in groups
        # Make one output CSV file for the whole event tree
        filename = "./outputs2/$init-groups.csv"
        output = ""
        for groupnum in 1:length(grouplist)
            # if init == "ISA"
            #     println("ISA group ", groupnum, " printing to CSV.")
            # end
            count = 0
            topseq = nothing
            # Find events that occur in all sequences of the group
            common_events = find_repeats(grouplist[groupnum])
            for seq in eachrow(grouplist[groupnum])
                count += 1
                initiator_BE = BEs[seq[:Init]]
                seqheader = string(", ", seq[:Seqnum],", ", initiator_BE.name,", ", initiator_BE.probability, ", ", seq[:Freq])
                # if init == "ISA"
                #     println(" Sequence ", seq[:Seqnum])
                # end
                if count == 1 # first sequence: print every event, even the
                              # shared ones
                    output *= string(groupnum)
                    seqheader *= string(", ", nrow(grouplist[groupnum]), " sequences, ", sum([s[:Freq] for s in eachrow(grouplist[groupnum])]))
                    output *= seqheader*"\n"
                    topseq = seq
                    for evt in seq[:Evts]
                        highlight = ""
                        if evt.probability <= 0.5
                            highlight *= "fail "
                            if !(evt.name in common_events)
                                highlight *= "extra "
                            end
                        end
                        if evt.name in common_events
                            highlight *= "repeat "
                        end
                        output *= string(highlight, ", , ", evt.name, ", ", evt.probability, ", , ", evt.description, "\n")
                    end
                else # all other sequences: print only events not shared by the
                    # whole group
                    # Check if this sequence is a non-minimal version of any
                    # of the previous sequences in this group
                    nms = nonminimal(seq, grouplist[groupnum], groupnum)
                    output *= seqheader*"\n"
                    for evt in seq[:Evts]
                        if !(evt.name in common_events) || evt.probability <= 0.5
                            highlight = ""
                            if evt.probability <= 0.5
                                highlight *= "fail "
                                if !(evt.name in common_events)
                                    highlight *= "extra "
                                end
                            end
                            if evt in nms[2]
                                highlight *= "nm "
                            end
                            output *= string(highlight, ", , ", evt.name, ", ", evt.probability, ", , ", evt.description, "\n")
                        end
                    end
                    output *= nms[1]
                end
            end # for seq in eachrow
            output *= "\n\n"
        end # for groupnum in grouplist
        # write the whole set of groups for the event tree to the file at once
        open(filename, "w") do file
            print(file, output)
        end
    end # for init, grouplist in groups
end # print_to_csv()

function main()
    # read all the sequence output from CSV files
    lines = readfiles()
    # separate it into sequences and the associated events
    BEs, inits, seqs = process_lines(lines)
    # separate them out by initiator
    seqs_by_init = group_by_init(seqs, inits)
    # Within each initiator, sort the sequences by key failure event.
    # Within each such group (of sequences that all share the same key failure)
    #  sum the sequence frequencies. Then sort the groups by total frequency
    #  to determine which groups are worth examining more closely.
    groups = group_by_failure(seqs_by_init)

    # This line creates the CSV output, which is only needed if
    # it hasn't been done before. Commented out to reduce runtime.
    #print_to_csv(groups, BEs)
    return groups, seqs_by_init, BEs
end

if ! @isdefined firstrun
    groups, seqs_by_init, BEs = main()
    firstrun = true
end

# Show the differences between two sequences
# First list is the events removed from A to get B
# Second list is the events added to A to get B
function seqdiff(seqA, seqB, verbose=false, failonly=false)
    diffs = (setdiff(seqA[:Evts], seqB[:Evts]), setdiff(seqB[:Evts], seqA[:Evts]))
    if failonly # remove any high-probability events from the list of
        # differences, since they're not key failures
        filter!(x->x.probability <= 0.94, diffs[1])
        filter!(x->x.probability <= 0.94, diffs[2])
        ## alternative: remove success events by name
        ## but leave in dependent failures
        # filter!(x->x.name[4] != 'A', diffs[1])
        # filter!(x->x.name[4] != 'A', diffs[2])

        ## a worse way of doing the same thing:
        # del1 = []
        # for f in diffs[1]
        #     if f.probability > 0.94
        #     # if f.name[4] == 'A'
        #         push!(del1, f)
        #
        #     end
        # end
        # del2 = []
        # for f in diffs[2]
        #     if f.probability > 0.94
        #     # if f.name[4] == 'A'
        #         push!(del2, f)
        #     end
        # end
        # for f in del1
        #     deleteat!(diffs[1],diffs[1] .== f)
        # end
        # for f in del2
        #     deleteat!(diffs[2],diffs[2] .== f)
        # end
    end
    len1 = length(diffs[1])
    len2 = length(diffs[2])
    if verbose
        if len1 == len2 == 0
            println("   Same.")
            return diffs
        else
            println("      Subtracted: ", "Added:")
        end
    end
    # format printing based on which list is longer
    for i in 1:max(len1,len2)
        n1 = "      "
        n2 = "      "
        if i <= len1
            d1i = diffs[1][i]
            n1 = d1i.name
        end
        if i <= len2
            d2i = diffs[2][i]
            n2 = d2i.name*"   "*d2i.description
        end
        if verbose
            println("      ", n1, "       ", n2)
        end
    end
    return diffs
end

#function combineseqs(seqs, sublists)


#function combine_seqs(grouplist, correspondences)

function compare_to_group(seq, group)
    group = eachrow(group)
    first = group[1]
    minlen = 1000
    comp = first
    for prev_seq in group[1:end]
        diffs = seqdiff(prev_seq, seq)
        len = length(diffs[1])+length(diffs[2])
        if len < minlen
            comp = prev_seq
            minlen = len
        end
    end
    println("Comparing ", seq[:Seqnum], " to ", comp[:Seqnum])
    seqdiff(comp, seq, true, true)
end

function compare_group(group, failonly=false)
    group = eachrow(group)
    first = group[1]
    for i in 2:length(group) #group[2:end]
        seq = group[i]
        minlen = 1000
        comp = first
        for prev_seq in group[1:i-1]
            diffs = seqdiff(prev_seq, seq)
            len = length(diffs[1])+length(diffs[2])
            if len < minlen
                comp = prev_seq
                minlen = len
            end
        end
        println(i, ": Comparing ", seq[:Seqnum], " to ", comp[:Seqnum])
        seqdiff(comp, seq, true, failonly)
    end
end

# Compares each sequence in group1 to its closest match in group2
function compare_groups(group1, group2, failonly=false)
    for seq in eachrow(group1)
        compare_to_group(seq, group2)
    end
    println("Group lengths: ", length(eachrow(group1)), " ", length(eachrow(group2)))
    sum1 = sum([s[:Freq] for s in eachrow(group1)])
    sum2 = sum([s[:Freq] for s in eachrow(group2)])
    println("Group frequencies: ", sum1, " ", sum2)
end

# Compares the first sequence in the group to each of the other sequences
function compare_to_first(group)
    group = eachrow(group)
    first = group[1]
    for i in 2:length(group) #group[2:end]
        seq = group[i]
        comp = first
        println("Comparing ", seq[:Seqnum], " to ", comp[:Seqnum])
        seqdiff(comp, seq, true, true)
    end
end

# never finished writing this function
# i've been using compare_group or compare_to_first instead
function describe_group(group)
    descrips = []
    group = eachrow(group)
    first = group[1]
    for i in 2:length(group) #group[2:end]
        seq = group[i]
        minlen = 1000
        comp = first
        for prev_seq in group[1:i-1]
            diffs = seqdiff(prev_seq, seq, false, true)
            len = length(diffs[1])+length(diffs[2])
            if len < minlen
                comp = prev_seq
                minlen = len
            end
        end
        diffs = seqdiff(comp, seq)
        println("\nSequence ", seq[:Seqnum], " is similar to ", comp[:Seqnum], " except that it:")
    end
end

# for use in endstuff
function sum_nonmin(group)
    tot = 0
    for seq in eachrow(group)
        if length(nonminimal(seq, group, 1)[1]) > 0
            tot += seq[:Freq]
        end
    end
    return tot
end

# ending / endstuff is a function to calculate the numbers that go at the end
# of the GPSA sequence comparison table for each event tree
function endstuff(name, numgroups, seqs_by_init, groups)
    println("Listed freq ", sum([g.Freq_sum for g in eachrow(by(seqs_by_init[name], :Failure, :Freq => sum))[1:numgroups]]))
    println("Nonmin ", sum([sum_nonmin(groups[name][i]) for i in 1:numgroups]))
end

ending = (n, ng) -> endstuff(n, ng, seqs_by_init, groups)

# This block is mostly a place to stick lines I may want to run in the REPL when
# I need specific info.
if firstrun # will always be true b/c i set it earlier, i think

    train10 = ["VETF10", "VEDA01", "TXTA09", "TFDA01", "TDTA07"]
    train30 = ["VETF30", "VEDA03", "TXTA01", "TFDA55", "TDTA01"]
    norev = ["V2DA01", "TXTA09", "HPQA01", "TPQA01", "RCQA01", "RHTA01"]
    rev20 = ["V2DC05", "TXTC11", "HPQE05", "TPQE05", "RCQE05", "RHTC03"]
    rev40 = ["V2DB02", "HPQB02", "TPQB02", "RCQB02"]
    for i in 1:length(train10)
        @assert BEs[train10[i]].probability == BEs[train30[i]].probability
    end

    group1 = eachrow(groups["SLBIC"][1])

    #seqdiff(eachrow(groups["LOVH"][2])[2], eachrow(groups["LOVH"][2])[7], true);

    #sum([sum_nonmin(groups["LOVH"][i]) for i in 1:4])

    #compare_group(groups["LOVN"][3])

    #BEs["TLF1"].description

    slbic_msrvo_1 = append!(copy(groups["SLBIC"][1]), copy(groups["MSRVO"][1]))

    # for seq in eachrow(groups["LOCV"][1])
    #    compare_to_group(seq, groups["LOVH"][1])
    # end

    plomfw_groupcdfs = by(seqs_by_init["PLOMFW"], :Failure, :Freq => sum);
    sum([g.Freq_sum for g in eachrow(plomfw_groupcdfs)[1:3]])

    sum([g.Freq_sum for g in eachrow(by(seqs_by_init["C3MSIV"], :Failure, :Freq => sum))[1:3]])

    # for seq in eachrow(groups["C3MSIV"][3])
    #     println(seq[:Seqnum])
    # end

    seqs_by_init["L400KV"][(seqs_by_init["L400KV"][!,:Seqnum] .== 9760),:][:Evts]
    # better version of that
    @where(seqs_by_init["L400KV"], :Seqnum .== 9760)[:Evts][1]

    #print(replace("=IF(INDEX(ISA!C:C,MATCH(A6,ISA!A:A,0))>0,INDEX(ISA!C:C,MATCH(A6,ISA!A:A,0)),\"\")","ISA"=>"LOWKI"))
end
