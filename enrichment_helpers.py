import numpy as np
import scipy.stats as stats
import multiprocessing as mp
import random
import ngramNets
from tqdm import tqdm

# Defining a function to calculate cohen's D since I do not want to make the assumption that there is equal variance between the random and subsampled distributions. 
def cohen_d(a,b):
    '''
    Calculates the Cohen's d effect size of 2 lists of values, where the second is considered the "control" group.
    '''

    
    if all(np.isnan(a)) and all(np.isnan(b)): # For certain cases there might be all nans so just check that first
        d = np.nan
    else:
        # Means
        x1 = np.nanmean(a)
        x2 = np.nanmean(b)

        # Sizes
        n1 = len(a)
        n2 = len(b)

        # Defining the pooled standard deviation (Note: need to account for the definition by cohen that uses N-1 for the variance, which numpy does not use)
        s = np.sqrt(((n1-1)*np.nanstd(a,ddof=1)**2 + (n2 -1)*np.nanstd(b, ddof=1)**2)/(n1+n2-2))
        if s > 0:
            d = (x1-x2)/s
        elif np.isnan(s):
            d = np.nan
        else:
            d = 0
        

    return d

def hypergeom_prune_ns(degnn, sweep):
    
    up_hyper_vals = degnn.ngram_DEG_hypergeom('Up')
    dn_hyper_vals = degnn.ngram_DEG_hypergeom('Down')
    ns_sweep = []
    
    for i in sweep:
        up_check = [k for k,v in up_hyper_vals.items() if v <= i]
        dn_check = [k for k,v in dn_hyper_vals.items() if v <= i]

        # Setting up the networks for determining the network separation
        up_net = degnn.G.subgraph(up_check)
        dn_net = degnn.G.subgraph(dn_check)
        if len(up_net.nodes()) > 0 and len(dn_net.nodes())> 0:
            ns = ngramNets.network_separation(up_net, dn_net, degnn.ref_data)
        else:
            ns = np.nan    
        ns_sweep.append(ns)

    return ns_sweep

# Below is a handful of helper functions to try and do network separation on a number of iterations to create a null distribution of DEG network separation values
def get_random_count_dist(count_weights):
    rand_counts = []
    for i in count_weights[1:11]:
        i_min = np.floor(i*.85)
        i_max = np.ceil(i*1.15)
        if i_min == 0 and i_max == 0:
            r = 0
        else:
            r = random.randrange(i_min, i_max)
        rand_counts.append(r)
    if count_weights[11] > 0:
        rand_counts.append(random.randrange(0,count_weights[11]))
    else:
        rand_counts.append(0)
    return rand_counts

def get_random_id_indices(arch_len_array, rand_count_list):
    rand_4_analysis = []
    for i in range(1,11):
        cands = arch_len_array[arch_len_array == i]
        rand_ids = random.sample(list(cands.index), k=rand_count_list[i-1])
        rand_4_analysis += rand_ids

    # For the >10 n-grams now getting candidates
    cands = arch_len_array[arch_len_array > 10]
    rand_ids = random.sample(list(cands.index), k=rand_count_list[10])
    rand_4_analysis += rand_ids

    return rand_4_analysis

def designate_rand_DEGs(rand_indices, ref_df, up_fraction):
    # Now double-checking the distribution
    random_prots = ref_df.loc[rand_indices]

    # Now designating the random genes as either up or down based on the proportion of DEGs that were for both
    rand_up_num = round(up_fraction*len(random_prots))

    # Alternative to the random sampling which has issues with reproducibility when called multiple times
    rand_prot_set = sorted(random_prots['UniProt ID'].tolist())
    #random.shuffle(rand_prot_set)
    #rand_up_chosen = rand_prot_set[0:rand_up_num]
    #rand_down_chosen = rand_prot_set[rand_up_num:]

    rand_up_chosen = sorted(random.sample(rand_prot_set, k=rand_up_num))
    rand_down_chosen = sorted(set(rand_prot_set).difference(rand_up_chosen))
    rand_DEGs = {'Up':rand_up_chosen,'Down':rand_down_chosen}
    return rand_DEGs

# Defining a helper function for grabbing the n-grams of interest for a hypergeometric p-value of interest.
def get_hyper_ngrams(degnn, cutoff):
    up_hyper = degnn.ngram_DEG_hypergeom('Up')
    dn_hyper = degnn.ngram_DEG_hypergeom('Down')

    ngrams_2_check = []
    for i in [up_hyper, dn_hyper]:
        k = []
        for node in i:
            if i[node] < cutoff:
                k.append(node)
        ngrams_2_check.append(k)
    
    return ngrams_2_check

def get_hyper_prots(degnn, ngram_lists):

    hyper_prots = []
    for i,cond_degs in zip(ngram_lists,[degnn.up_DEGs, degnn.down_DEGs]):
        k = set()
        for ngram in i:
            j = degnn.interpro2uniprot[ngram]
            inter = set(j).intersection(cond_degs)
            k.update(inter)
        hyper_prots.append(list(k))

    return hyper_prots


def get_random_net_sep_metrics(normalized_arch_weights, complete_arch_dist, degnn,deg_ratio,pvals= np.logspace(0,-10, num=21),return_dist = False):

    rcd = get_random_count_dist(normalized_arch_weights)
    rii = get_random_id_indices(complete_arch_dist,rcd)
    r_prots = designate_rand_DEGs(rii, degnn.ref,deg_ratio)
    
    # This is for debugging as necessary
    #print('Up Random:')
    #print(r_prots['Up'][0:10])
    #print('Down Random:')
    #print(r_prots['Down'][0:10])

    degnn.set_DEG_ngrams(r_prots['Up'],r_prots['Down'], verbose=False)
    ns = degnn.DEG_network_sep()
    ns_vals = hypergeom_prune_ns(degnn, pvals)
    ns_iqr = stats.iqr(ns_vals, nan_policy='omit')
    if return_dist:
            o = (ns,ns_iqr,ns_vals)
    else:
        o = (ns,ns_iqr)

    return o

def get_subsample_net_sep_metrics(degnn,orig_DEGs, rand_prot_nums,pvals= np.logspace(0,-10, num=21),return_dist = False):
    
    # Getting a random sample of the up DEGs (not accounting for the n-gram length)
    orig_up = orig_DEGs[0]
    orig_dn = orig_DEGs[1]

    rand_up = random.sample(orig_up, k=rand_prot_nums[0])
    rand_dn = random.sample(orig_dn, k=rand_prot_nums[1])

    # This is for debugging as necessary
    #with open('debug.log', 'a') as f:
    #    f.write('\nUp Subsample:')
    #    for p in rand_up[0:10]:
    #        f.write(f"{p}\t")
    #print('Down Subsample:')
    #print(rand_dn[0:10])

    degnn.set_DEG_ngrams(rand_up,rand_dn, verbose=False)
    ns = degnn.DEG_network_sep()
    ns_vals = hypergeom_prune_ns(degnn,pvals)
    ns_iqr = stats.iqr(ns_vals, nan_policy='omit')

    if return_dist:
        o = (ns,ns_iqr,ns_vals)
    else:
        o = (ns,ns_iqr)

    return o

def individual_trial_calc(degnn, arch_weights, comp_arch_dist, ratio, sweep, originals,dist_flag = False, seed =123):
    random.seed(seed)
    rand_data = get_random_net_sep_metrics(normalized_arch_weights=arch_weights,
                                          complete_arch_dist=comp_arch_dist,
                                          degnn=degnn,
                                          deg_ratio=ratio,
                                          pvals=sweep,
                                          return_dist=dist_flag)
    # Since the random function above generates the sizes will pass this on.
    rand_szs = (len(degnn.up_DEGs), len(degnn.down_DEGs))
    subsample_data = get_subsample_net_sep_metrics(degnn=degnn,
                                orig_DEGs=originals,
                                rand_prot_nums=rand_szs,
                                pvals=sweep,
                                return_dist=dist_flag)
    
    return list(rand_data + subsample_data)
    
def calculate_separation_stability(degnn, num_trials = 50, pval_sweep = np.logspace(0,-10,21), return_distributions = False,processes = 1, verbose=True, progress_bar = False):
    orig_up = degnn.up_DEGs
    orig_dn = degnn.down_DEGs
    original_DEGs = list(orig_up) + list(orig_dn)
    deg_info = degnn.retrieve_protein_info(original_DEGs)
    arch_lens = deg_info['Interpro Domain Architecture IDs'].apply(lambda x: len(x.split('|')))
    max_length = max([max(arch_lens),11])
    deg_arch_dist = np.histogram(arch_lens, bins=range(max_length+1))
    complete_arch_dist = degnn.ref['Interpro Domain Architecture IDs'].apply(lambda x: len(x.split('|')))
    weight_list = deg_arch_dist[0][0:11]
    weight_list = np.append(weight_list,np.sum(deg_arch_dist[0][11:]))

    # Now normalizing the list 70% of the total number of DEGs available (depends on relative fractions. Don't want to bother making this completely correct)
    weight_list_n = np.round(weight_list/sum(weight_list)*(sum(weight_list)*0.7))
    up_frac = len(degnn.up_DEGs)/len(deg_info)

    # For reproducibility (mostly only for when multiprocessing is enacted) creating a seed list that will be passed to the individual trials function
    seedlist = random.sample(range(50*num_trials), num_trials)
    if processes == 1:
        if progress_bar:
            pbar = tqdm(total=num_trials)

        # Intializing some of the variables of the for loop
        rand_ns = np.zeros(num_trials)
        subsample_ns = np.zeros(num_trials)
        rand_iqr = np.zeros(num_trials)
        subsample_iqr = np.zeros(num_trials)
        
        if return_distributions: 
            rand_ns_dists = []
            subsample_ns_dists = []

        for i in range(num_trials):
            a = individual_trial_calc(degnn, weight_list_n,complete_arch_dist,up_frac,pval_sweep,[orig_up,orig_dn],dist_flag=return_distributions, seed=seedlist[i])
            # Unpacking the values into the datastructures above
            rand_ns[i] = a[0]
            rand_iqr[i] = a[1]
            if return_distributions:
                rand_ns_dists.append(a[2])
                subsample_ns[i] = a[3]
                subsample_iqr[i] = a[4]
                subsample_ns_dists.append(a[5])
            else:
                subsample_ns[i] = a[2]
                subsample_iqr[i] = a[3]

            if progress_bar:
                pbar.update(1)


        if progress_bar:
            pbar.close()
    else:
        if verbose:
            print('Will do multiprocessing')

        # Setting up the arguments to be passed to multiple processes
        args = [(degnn, weight_list_n,complete_arch_dist,up_frac,pval_sweep,[orig_up,orig_dn],return_distributions, seedlist[i]) for i in range(num_trials)]
       
        if __name__ == 'hypergeom_helpers':
            pool = mp.Pool(processes=processes)
            with pool as p:
                a = p.starmap(individual_trial_calc, args,chunksize=5)
            pool.close() # In case

        # Unpacking the values into the datastructures above
        rand_ns = [x[0] for x in a]
        rand_iqr = [x[1] for x in a]
        if return_distributions:
            rand_ns_dists = [x[2] for x in a]
            subsample_ns = [x[3] for x in a]
            subsample_iqr = [x[4] for x in a]
            subsample_ns_dists= [x[5] for x in a]
        else:
            subsample_ns = [x[2] for x in a]
            subsample_iqr = [x[3] for x in a]

    

    # Now returning the desired outputs

    if return_distributions:
        o = [rand_ns,rand_iqr,rand_ns_dists,subsample_ns,subsample_iqr,subsample_ns_dists]
    else:
        o = [rand_ns,rand_iqr,subsample_ns,subsample_iqr]

    return o

def retrieve_fpr_checks(degnn,num_DEGs,fpr_trials = 50, num_internal_trials = 50, deg_ratios = 0.6, processes = 1, progress_bar = False, seed =123):
    
    # Setting up the random DEGs
    # Here setting up a random seed list that ensures reproducibility across experimental runs
    rand_DEGs = degnn.retrieve_random_ids(num=num_DEGs, iters=fpr_trials,seed = seed)

    internal_fpr = []
    if progress_bar:
        pbar = tqdm(total=fpr_trials)

    for i in range(fpr_trials):
        cur_DEGs = next(rand_DEGs)
        cur_DEGs = sorted(cur_DEGs) # In case
        orig_up = random.sample(cur_DEGs, k=round(len(cur_DEGs)*deg_ratios))
        orig_dn = sorted(set(cur_DEGs).difference(orig_up))
        degnn.set_DEG_ngrams(up_DEGs=orig_up, down_DEGs=orig_dn, verbose=False)
        temp = calculate_separation_stability(degnn,num_trials=num_internal_trials,processes=processes,verbose = False)
        rand_iqr = [x for x in temp[1]]
        subsample_iqr = [x for x in temp[3]]
        rand_full_ns = [x for x in temp[0]]
        actual_full_ns = [x for x in temp[2]]

        # Now getting some of the stats and results
        temp_iqr_res = stats.mannwhitneyu(rand_iqr,subsample_iqr)
        temp_ns_res = stats.mannwhitneyu(rand_full_ns, actual_full_ns)

        internal_fpr.append((temp_ns_res[1], temp_iqr_res[1]))

        if progress_bar:
            pbar.update()

    if progress_bar:
        pbar.close()


    return internal_fpr

def calculate_fpr(actual_res, random_res):

    ns_res = actual_res[0]
    iqr_res = actual_res[1]
    ns_rand_res = [x[0] for x in random_res]
    iqr_rand_res = [x[1] for x in random_res]

    # Now calculating the fpr values
    x = [i <= ns_res for i in ns_rand_res]
    ns_fpr = sum(x)/len(x)
    x = [i <= iqr_res for i in iqr_rand_res]
    iqr_fpr = sum(x)/len(x)

    return (ns_fpr, iqr_fpr)