network = struct([]);
for i = 1:length(graph)
    network(i).node_one = graph(i).node_one;
    network(i).node_two = graph(i).node_two;
    network(i).length = graph(i).length;
end

nodeList_cmpl = []; clear i;
for j = 1:length(graph)
    n1 = graph(j).node_one;
    n2 = graph(j).node_two; 
    cmpl1 = n1(1) + (n1(2)*i);
    cmpl2 = n2(1) + (n2(2)*i);
    if sum(nodeList_cmpl == cmpl1) == 0
        nodeList_cmpl(end+1) = cmpl1;
    end
    if sum(nodeList_cmpl == cmpl2) == 0
        nodeList_cmpl(end+1) = cmpl2;
    end
end

nodeList = [];
for j = 1:length(nodeList_cmpl)
    nodeList(j,1) = real(nodeList_cmpl(j));
    nodeList(j,2) = imag(nodeList_cmpl(j));
end




% Initialize a degree counter
numNodes = size(nodeList, 1);
degreeCounter = zeros(numNodes, 1);

% Count the degrees
for i = 1:length(graph)
    node1 = find(ismember(nodeList, graph(i).node_one, 'rows'));
    node2 = find(ismember(nodeList, graph(i).node_two, 'rows'));
    degreeCounter(node1) = degreeCounter(node1) + 1;
    degreeCounter(node2) = degreeCounter(node2) + 1;
end

% Identify nodes with degree 1
nodesToRemove = find(degreeCounter == 1);

% Remove nodes and corresponding edges
newGraph = graph;
for i = length(graph):-1:1
    node1 = find(ismember(nodeList, graph(i).node_one, 'rows'));
    node2 = find(ismember(nodeList, graph(i).node_two, 'rows'));
    if ismember(node1, nodesToRemove) || ismember(node2, nodesToRemove)
        newGraph(i) = [];
    end
end

% Update nodeList if necessary
nodeList(nodesToRemove, :) = [];


source = 'Data/sample three/sample three';

writetable(struct2table(network),[source '_network_edges.csv'])
dlmwrite([source '_network_nodes.txt'], nodeList, 'delimiter','\t')
save([source '_graph.mat'],'graph')
