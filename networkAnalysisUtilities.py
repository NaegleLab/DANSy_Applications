import networkx as nx
import numpy as np
import itertools
import random
import networkAnalysisUtilities

def bridging_centrality(G):
    '''
    Calcualtes the bridging centrality of every node within the network Graph G. Bridging centrality is defined in Hwang et al. 2007

    Parameters:
    -----------
        - G: NetworkX Graph
            Graph containing the nodes of interest
    
    Returns:
    --------
        - bridge_cent: dict
            Dict containing the computed bridging centrality value for each node
            
    '''
    # Getting the dict of all betweenness centrality values
    bridge_cent = {}
    bc = nx.betweenness_centrality(G)
    for node in nx.nodes(G):
        bc_node = bc[node]
        degree_node = nx.degree(G, node)
        neighbor_degrees = dict(nx.degree(G, nx.neighbors(G, node)))
        n_degrees = [v for v in neighbor_degrees.values()]
        n_degrees = np.array(n_degrees)
        bridge_cent[node] = (1/degree_node)/(sum(1/n_degrees))*bc_node

    return bridge_cent

def flower_module_detection(G):
    '''
    There are specific connected components that look like flowers when it consists of >5 nodes. In these modules there is a central hub node connected to several terminal nodes. Here using a metric similar to the connectivity coefficient will try to identify these modules.
    
    To metric used is taking the node with the highest degree and if all the neighbors have a degree of one designating it as a flower module.

    Parameters:
    -----------
        - G: NetworkX Graph

    Returns:
    --------
        - flower_modules: dict
            dict containing the connected component information
    
    '''

    flower_modules = {}
    for c in nx.connected_components(G):
        d = dict(nx.degree(G, c))
        d_list = [v for v in d.values()]
        d_array = np.array(d_list)

        # Small feature of these modules is that the degree of the hub node is equal to exactly half of the sum of the entire connected component
        if sum(d_array)/2 == max(d_array):
            # Finding the hub node
            for k,v in d.items():
                if v == max(d_array):
                    hub = k

            # Putting the full connected component and hub to be outputted with the hub being the identifier
            flower_modules[hub] = c

    return flower_modules

def wheel_module_detection(G):
    '''
    Finds and returns the connected components that are hub and spoke topologies (flowers as I call them), where each spoke is fully connected simliar to a wheel-like topology.

    Parameters:
    -----------
        - G: NetworkX Graph

    Returns:
    --------
        - wheel_modules: dict
            dict containing the connected component information with the hub node as the key
    '''

    wheel_modules = {}
    for c in nx.connected_components(G):
        d = dict(nx.degree(G, c))
        d_list = [v for v in d.values()]
        d_array = np.array(d_list)

        # Small feature of these modules is that the degree of the hub node is equal to exactly a quarter of the sum of the entire connected component. However, there a square topology between 4 nodes will also be identified. To differentiate just need to check if the minimum and maximum node values are the same.
        if sum(d_array)/4 == max(d_array) and max(d_array) != min(d_array):
            # Finding the hub node
            for k,v in d.items():
                if v == max(d_array):
                    hub = k

            # Putting the full connected component and hub to be outputted with the hub being the identifier
            wheel_modules[hub] = c

    return wheel_modules

def ring_module_detection(G):
    '''
    Finds and returns the connected components that ring like structures where all nodes have a degree of 2 and equal betweenness centrality values for each edge.

    Parameters:
    -----------
        - G: NetworkX Graph

    Returns:
    --------
        - ring_modules: dict
            dict containing the connected component information
    '''

    ring_modules = {}
    for c in nx.connected_components(G):
        d = dict(nx.degree(G, c))
        d_list = [v for v in d.values()]
        all_equal = len(set(d_list)) == 1

        if all_equal and d_list[0] == 2:
            nodes_c = [n for n in d.keys()]
            edge_cent = nx.edge_betweenness_centrality_subset(nx.subgraph(G,nodes_c), sources=nodes_c, targets=nodes_c)
            edge_cent_equal = len(set([v for v in edge_cent.values()])) == 1
            if edge_cent_equal:
                # To name the community using the node with the longest n-gram as it likely contains the most information
                n_length = {x:len(x.split('|')) for x in nodes_c}
                hub = max(n_length, key=n_length.get)
                ring_modules[hub] = c

    return ring_modules

def contract_simple_cc_modules(G):
    '''
    Contracts a network by collapsing connected components that show flower-like or wheel-like topologies.

    Parameters:
    -----------
        - G: NewtworkX Graph

    Returns:
    --------
        - G_contracted: NetworkX Graph
            The new contracted network with modules as isolates
        - cc_communities: dict
            dict containing the isolate node name and the member nodes that were collapsed into it
    '''

    G_contracted = G.copy()
    cc_communities = {}

    flowers = flower_module_detection(G_contracted)
    wheels = wheel_module_detection(G_contracted)
    rings = ring_module_detection(G_contracted)

    cc_communities = {**flowers, **wheels,**rings}

    for module, members in cc_communities.items():
        for node in members:
            if node != module:
                nx.contracted_nodes(G_contracted, module, node, False,False)

    return G_contracted, cc_communities


def contract_network_via_GN_communities(G):
    '''
    Takes a networkX graph and finds communities as defined through the Girvan-Newman algorithm. It then takes these communities and contracts the nodes into the node with the highest degree and any nodes that help connect communities.

    Parameters:
    -----------
        - G: NetworkX Graph

    Returns:
    --------
        - G_contracted: NetworkX Graph
            The new network that has been contracted. This is a copy of the original.
        - community_OI: list
            The specific communities identified in the network
    '''

    # Finding the communities via the GN algorithm using the betweenness centrality. Note: Future work may allow this to vary in which centrality measure will be used.
    communities = nx.community.girvan_newman(G)
    modularity_check = 0
    
    # Getting the communities that have the highest modularity as a starting cutoff. (Will have to look into other potential parameters to help select this.)
    if G.number_of_edges() == 0:
        community_OI = list(nx.connected_components(G))
    else:
        for comm in itertools.takewhile(lambda c: nx.community.modularity(G,c) > modularity_check, communities):
            modularity_check = nx.community.modularity(G,comm)
            last_community = comm
        community_OI = last_community

    G_contracted = G.copy()
    node_community_degrees = retrieve_node_community_degree(G_contracted, community_OI)
    for c in community_OI:
        deg_comm = dict(nx.degree(G_contracted,c))
        hub_node = max(deg_comm, key=deg_comm.get) # This can create potential issues later if more than one node has the highest degree
        for node in c:
            if node is not hub_node and node_community_degrees[node]['K_out'] == 0:
                nx.contracted_nodes(G_contracted, hub_node, node, False, False)
    
    return G_contracted, community_OI
            

def retrieve_node_community_degree(G, community):
    ''''
    Calculates both the intra-community and inter-community degrees for nodes within individual communities.

    Parameters:
    -----------
        - G: NetworkX Graph
            Network that the nodes and communities originate from
        - c: list
            Contains the individual communities as individual elements and their corresponding nodes

    Returns:
    --------
        - node_comm_deg: dict
            dict with the intra- and inter-community degrees for each node
    '''

    # Define the community the node is a part of for quick assessment of whether edges are within or outside a community
    node_community = {}
    for i, c in enumerate(community):
        for node in c:
            node_community[node] = i

    # Intialize the dict with 0 degrees for each node, K_in is within the community, K_out is outside the community
    node_comm_deg = {}
    for node in nx.nodes(G):
        node_comm_deg[node] = {'K_in':0, 'K_out':0}

    for e in nx.edges(G):
        if node_community[e[0]] == node_community[e[1]]:
            node_comm_deg[e[0]]['K_in'] += 1
            node_comm_deg[e[1]]['K_in'] += 1
        else:
            node_comm_deg[e[0]]['K_out'] += 1
            node_comm_deg[e[1]]['K_out'] += 1
    
    return node_comm_deg

def find_community_strength(G, communities):
    '''
    Determines whether a community is considered a weak or strong community based on the inter- and intra- community degrees as outlined in Radicchi, et al PNAS 2004 (https://doi.org/10.1073/pnas.040005410)

    Parameters:
    -----------
        - G: networkX Graph
            Network that the communities originated from
        - communities: list
            Contains a list of dicts for each community with a potential name and their member nodes

    Returns:
    --------
        - community_strength: dict
            dict containing communities, their strenght, and a potential name based on the hub node
    '''
    all_comms = [c['Community'] for c in communities]
    node_degrees = retrieve_node_community_degree(G, all_comms)
    community_strength = {}
    current_c = 0
    for c in communities:
        deg_comm = dict(nx.degree(G,c['Community']))
        hub_node = c['Potential Name']
        weak_community_flag = 0
        for node in c['Community']:
            # This checks if a node is more strongly connected to the community versus outside it.
            if node_degrees[node]['K_in'] <= node_degrees[node]['K_out']:
                weak_community_flag = 1
                break
        
        if weak_community_flag:
            inter_degree_sum = 0
            intra_degree_sum = 0
            for node in c['Community']:
                inter_degree_sum += node_degrees[node]['K_out']
                intra_degree_sum += node_degrees[node]['K_in']
            
            if intra_degree_sum > inter_degree_sum:
                community_status = 'Weak'
            elif len(c['Community']) == 1 and inter_degree_sum == 0:
                community_status = 'Isolate'
            else:
                community_status = 'Uh oh'
        else:
            community_status = 'Strong'
        
        community_strength[current_c] = {'Community':c['Community'],
                              'Potential Name':hub_node,
                              'Strength':community_status}
        current_c += 1
    
    return community_strength

def contract_network_via_louvian_communities(G, seed):
    '''
    Takes a networkX graph and finds communities as defined through the Louvian algorithm. It then takes these communities and contracts the nodes into the node with the highest degree and any nodes that help connect communities.

    Parameters:
    -----------
        - G: NetworkX Graph

    Returns:
    --------
        - G_contracted: NetworkX Graph
            The new network that has been contracted. This is a copy of the original.
        - community_list: list
            The specific communities identified in the network
    '''
    if not seed:
        seed = None
    
    # Finding the communities via the GN algorithm using the betweenness centrality. Note: Future work may allow this to vary in which centrality measure will be used.
    community_OI = nx.community.louvain_communities(G, weight = None, seed=seed)
    community_list = []
    G_contracted = G.copy()
    node_community_degrees = retrieve_node_community_degree(G, community_OI)
    
    for c in community_OI:
        deg_comm = dict(nx.degree(G,c))
        hub_node = max(deg_comm, key=deg_comm.get) # This can create potential issues later if more than one node has the highest degree
        for node in c:
            if node is not hub_node and node_community_degrees[node]['K_out'] == 0:
                nx.contracted_nodes(G_contracted, hub_node, node, False, False)
        community_list.append({'Potential Name':hub_node,
                               'Community':c})
    
    return G_contracted, community_list

### Functions used for determing community robustness and overlap over time 


def create_iteration_seeds(k = 100, seed = []):
    '''
    Generates a list of random integer seeds for the k iterations.

    Parameters:
    -----------
        - k: int
            Number of iterations to generate seeds for
        - seed: int (Optional)
            Seed for the RNG. If not provided will use the datestamp.

    Returns:
    --------
        - x: list
            list of randomly distributed and unique integers
    '''
    
    if seed:
        random.seed(seed)
    else:
        random.seed()

    # Ensuring enough unique random numbers are generated
    if k > 1000:
        max_k = k*2
    else:
        max_k  = 1000

    x = random.sample(range(max_k), k)

    return x

def generate_community_reference(G, seed_list):
    '''
    Generates a dictionary containing a reference set of community detection iterations of the Louvian algorithm based on the inputted random seeds.

    Parameters:
    -----------
        - G: networkX Graph
            - Graph that communities will be detected from
        - seed_list: list
            - list of integers that will be used for the Louvian algorithm

    Returns:
    --------
        - comm_checks: dict
            - dictionary of each iteration consisting of the community detection details including the community number, seed, and community details as values
    '''
    
    comm_checks = {}
    # Ensuring each iteration has it's own entry
    cnt = 0
    for seed in seed_list:
        _, rand_comms = networkAnalysisUtilities.contract_network_via_louvian_communities(G, seed = seed)
        comm_checks[cnt] = {'Seed':seed,
                            'Communities': rand_comms,
                            'Number Communities':len(rand_comms)}
        cnt += 1
    return comm_checks

def extract_comm_iter_info(community_iters,k):
    '''
    Extract community iteration information about what are core and peripheral members of each community, the number of instances it occurs, and whether it was consistently identified.

    Parameters:
    -----------
        - community_iters: dict
            dictionary containing the details of communities found in each iteration of the community detection algorithm

    Returns:
    --------
        - community_ii: dict
            dictionary containing summarizing details of any community detected across all iterations
    '''

    community_ii = {}
    
    for community_iteration in community_iters:
        
        comms = community_iters[community_iteration]['Communities']
        
        # Going through each community in the specific iteration
        for c in comms:
            comm_name = c['Potential Name']
            comm_members = c['Community']
            
            if comm_name in community_ii:
               
                # Putting the different iterations (the plan is for these to be compared to each other for similarity)
                community_ii[comm_name]['Iterations'].append(comm_members)
                community_ii[comm_name]['Instances'] += 1
                
                # Checking core (consistently identified) n-grams of the communities
                old_core = community_ii[comm_name]['Core']
                core = old_core.intersection(comm_members)
                
                if old_core != core:
                    community_ii[comm_name]['Core'] = core
                    community_ii[comm_name]['Changed Flag'] = 1

                    # If there is a change keeping track of those which were dropped
                    core_peripheral = old_core.difference(core)
                else:
                    core_peripheral = set() #empty set for later
                
                # Finding the non-core n-grams and designating them as peripheral
                peripheral = comm_members.difference(core)

                if 'Peripheral' in community_ii[comm_name]:
                    
                    # Ensuring the peripheral n-grams are retained for downstream analysis
                    old_peripheral = community_ii[comm_name]['Peripheral']
                    if peripheral != old_peripheral or core_peripheral:
                        peripheral = peripheral.union(old_peripheral).union(core_peripheral)
                        community_ii[comm_name]['Peripheral'] = peripheral
                
                # If the core changed or the community that was checked has more than the core create the peripheral n-grams
                elif peripheral or core_peripheral:
                    community_ii[comm_name]['Peripheral'] = peripheral.union(core_peripheral)

            # If the first time a community is detected  
            else:
                community_ii[comm_name]= {'Iterations':[comm_members],
                                            'Core':comm_members,
                                            'Changed Flag':0,
                                            'Instances':1}
            
    # Now that all the information for each community has been found will be going through and ensuring the communities do refer to different ones.
    community_ii = remove_duplicate_communities(community_ii,k)

    return community_ii

def remove_duplicate_communities(community_info,k):
    '''
    Determines whether any communities are actually the same community with different names. This becomes necessary due to earlier decision to use the name of the node with the highest degree for the community name.

    Parameters:
    -----------
        - community_info: dict
            - dictionary of communities and their iterations and instance number

    Returns:
    --------
        - community_info: dict
            - modified version of the input dictionary where the communities with different names are collapsed into each other and modified accordingly.
    
    '''

    # Grab any potential communities that could be duplicates. These should only occur if the community has not had its core members change at all. (Future me will have to check how true this is and change it accordingly.)
    potential_dups = {c:v['Core'] for c,v in community_info.items() if v['Instances'] < k and v['Changed Flag'] == 0}
    dups = []

    # Run through all the potential duplicates and check if their core members match those of another community.
    for c, core in potential_dups.items():
        for c_2, core_2 in potential_dups.items():
            if c != c_2: # Don't check self
                core_comp = set(core).difference(core_2)
                core_comp_i = set(core_2).difference(core) # Check the inverse.
                if (len(core_comp) == 0) and (len(core_comp_i) == 0): # Both have to be exactly zero. If want to check if one is a subcommunity then can change to an or statement
                    dups.append([c,c_2]) # List of pairs

    # Run through each pair of communities and combine them followed by removing the other one.
    for pair in dups:
        # Check that both communities are still present
        c1 = pair[0]
        c2 = pair[1]
        
        if (c1 in community_info) and (c2 in community_info):
            #print('Combining %s and %s'%(c1, c2)) # Remove later
            # Determine which name is shorter as it likely is more informative of potential connections
            c1_name_len = len(c1.split('|'))
            c2_name_len = len(c2.split('|'))
            c1_num = community_info[c1]['Instances']
            c2_num = community_info[c2]['Instances']
            if c1_name_len > c2_name_len:
                community_info[c1]['Instances'] += c2_num
                community_info[c1]['Iterations'].append(community_info[c2]['Iterations'])
                del community_info[c2]
            else:
                community_info[c2]['Instances'] += c1_num 
                community_info[c2]['Iterations'].append(community_info[c1]['Iterations'])
                del community_info[c1]

    return community_info

def simplify_network(G):
    '''
    Contracts a network to remove isolates and the simple connected component topologies. This is meant to only be used for reducing computational load of any community detection algorithms.

    Parameters:
    -----------
        - G: networkX Graph
            - Graph to be simplified
    
    Returns:
    --------
        - G_contracted: networkX Graph
            - Contracted version of the inputted network
    '''

    # Make copy and perform all operations on this in case the original is wanting to be kept.
    G_contracted = G.copy()

    # Finding and removing isolates
    isos = list(nx.isolates(G_contracted))
    G_contracted.remove_nodes_from(isos)
    G_contracted, _ = contract_simple_cc_modules(G_contracted)
    simple_mods = list(nx.isolates(G_contracted))
    G_contracted.remove_nodes_from(simple_mods)

    return G_contracted

def compare_soft_cluster_freq(ref_comm, comp_comm, thres = 0, q = 0.99):
    '''
    Calculates the difference in the frequency distribution of node membership to different soft cluster communities between the two inputs. A threshold to designate differences as potential noise or not can be provided. If not provided then it is calculated based on the input quantile (or 99% if omitted).

    Parameters:
    -----------
        - ref_comm: dict
            - dictionary containing the soft clustering results for the reference dataset
        - comp_comm: dict
            - dictionary containing the soft clustering results for the comparative dataset
        - thres: float (Optional)
            - threshold cutoff for designating frequency differences as noise
        - q: float (Optional)
            - quantile fraction to determine the threshold

    Returns:
    --------
        - diff_communities: list
            - list of the community names found to be different from the reference dataset
    '''


    # Extracting the differences in the frequencies for each node
    temp_comp = calc_soft_clust_membership_diff(comp_comm, ref_comm)
    
    # Threshold is either an input or determined via the distribution of all values
    if thres == 0:
        abs_diff = []
        for c_mem in temp_comp.values():
            for x in c_mem:       
                abs_diff.append(abs(c_mem[x]))
        thres = np.quantile(abs_diff, q)

    # Generating a community-centric colleciton of the frequency differences
    comm_freq_diff = {}
    for node, f_details in temp_comp.items():
        for c in f_details.keys():
            if c not in comm_freq_diff:
                comm_freq_diff[c] = [f_details[c]]
            else:
                comm_freq_diff[c].append(f_details[c])

    # Filtering out the communities to those that had some difference across all iterations.
    c_2_keep = [k for k,v in comm_freq_diff.items() if not any(v)]
    for c in c_2_keep:
        del comm_freq_diff[c]

    # Determining which nodes are noise
    

    noise_comms = []
    for node, f in comm_freq_diff.items():
        if all([abs(v) <= thres for v in f]):
            noise_comms.append(node)

    # Removing the frequency differences that were associated with noise only
    for c in noise_comms:
        del comm_freq_diff[c]

    # The names of communities that are different
    diff_communities = [k for k in comm_freq_diff.keys()]

    return diff_communities

def calc_soft_clust_membership_diff(comm_1, comm_2):
    '''
    Calculates the difference in frequencies from a reference dataset for a node's membership in different communities via soft clustering.

    Parameters:
    -----------
        - comm_1: dict
            - Frequency distributions for the comparative dataset
        - comm_2: dict
            - Frequency distributions for the reference dataset
    
    Returns:
    --------
        - temp_comp: dict
            - Differences in frequency distributions
    '''
    
    temp_comp = {}
    
    # For each node in the comparative dataset retrieve the community and frequency values
    for node in comm_1:
        
        # Temporary dictionary containing frequencies for the node in the comparative dataset
        tmp = {}
        for x in enumerate(comm_1[node]['Communities']):
            tmp[x[1]] = comm_1[node]['Frequency'][x[0]]
        
        # Temporary dicitonary containing frequencies for the node in the reference dataset
        tmp2 = {}
        if node in comm_2:
            for x in enumerate(comm_2[node]['Communities']):
                tmp2[x[1]] = comm_2[node]['Frequency'][x[0]]

        # Determining the differences and ensuring all communities detected across both datasets are accounted for
        temp_comp[node] = {}
        for x in list(set(tmp.keys()).union(tmp2.keys())):
            if (x in tmp2) and (x in tmp):
                temp_comp[node][x] = tmp[x] - tmp2[x]
            elif x in tmp:
                temp_comp[node][x] = tmp[x]
            elif x in tmp2:
                temp_comp[node][x] = -tmp2[x]

    return temp_comp

def get_freq_diff_noise_thres(ref_comm_1, ref_comm_2, q = 0.99):
    '''
    Calculate and get the threshold for potential noise given the reference dataset of soft clustering results.

    Parameters:
    -----------
        - ref_comm_1: dict
            - Frequency distribution of one set of reference communities
        - ref_comm_2: dict
            - Frequency distribtuion of a second set of reference communities
        - q: float (Optional)
            - quantile that the threshold should be determined as
    
    Returns:
    --------
        - thres: float
            - Threshold for the frequency difference distributions
    '''

    ref_comp = calc_soft_clust_membership_diff(ref_comm_1, ref_comm_2)
    abs_diff = []
    for c_mem in ref_comp.values():
        for x in c_mem:       
            abs_diff.append(abs(c_mem[x]))
    thres = np.quantile(abs_diff, q)

    return thres

def get_node_soft_clust_res(comm_for_check, k):
    '''
    Returns the soft cluster results for each node within a network. (Note: this may become deprecated in the future by calculating frequency earlier in the process.)

    Parameters:
    -----------
        - comm_for_check: dict
            - dictionary containing all the community information across all iterations
        - k: int
            - number of iterations used for community detection

    Returns:
    --------
        - node_freq: dict
            - Frequency distribution of soft cluster results for each node
    '''

    node_freq = {}
    for comm, details in comm_for_check.items():
        insta = details['Instances']
        if 'Peripheral' in details:
            for node in details['Peripheral']:
                iter_counts = 0
                
                # Going through each iteration and determining how often the node occurs
                for iter in details['Iterations']:
                    if node in iter:
                        iter_counts += 1
                freq = (iter_counts/insta) * (insta/k)
                if node not in node_freq:
                    node_freq[node] = {'Communities':[comm],'Frequency':[freq]}
                else:
                    node_freq[node]['Communities'].append(comm)
                    node_freq[node]['Frequency'].append(freq)
        
        # Core nodes are found consistently so these are associated with how often the community is found instead
        for node in details['Core']:
            if node not in node_freq:
                node_freq[node] = {'Communities':[comm],'Frequency':[insta/k]} #Scaling based on how often the community is detected
            else:
                node_freq[node]['Communities'].append(comm)
                node_freq[node]['Frequency'].append(insta/k)
    
    return node_freq