 alphabet = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n'
'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
 p = [ .08167 .01492 .02782 .04253 .12702 .02228 .02015 .06094
.06966 .00153 .00772 .04025 .02406 .06749 .07507 .01929 .00095 .05987
.06327 .09056 .02758 .00978 .02361 .00150 .01974 .00074 ];
 alpha=my_epektash(alphabet);
 
 pr=my_prob(p);
 dict=mey_huffmandict(alpha,pr);
 sourcebhta=importdata('kwords.txt');
 sourcebhta=[sourcebhta{:}];
 sourceb=lower(sourcebhta);
 sigb={};
 for i=1 : length(sourceb)-1
 if ismember(lower(sourceb(i)),alphabet)
 sigb=[sigb strcat(sourceb(i),sourceb(i+1))];
 end
 end


 encob=myhuffmanenco(sigb,dict);
 dsigb = my_huffmandeco(encob,dict);
 entropia=0;
 for l=1 : 676
 entropia= entropia +pr(l)*log2(1/pr(l));
 end
 mesos={};

 for j=1 :676
 mesos(j) = dict(j,2);

 end
 length=cellfun('length',mesos);
 true_mesos=0;
 for k=1 : 676
 true_mesos= true_mesos + pr(1,k)*length(1,k)';
 end

 apodosh=entropia/true_mesos;
%ERWTHMA 5b 
 alphabet = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n'
'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
 sourcebhta=importdata('kwords.txt');
 sourcebhta=[sourcebhta{:}];
 sourceb=lower(sourcebhta);
 sigb={};
 for i=1 : length(sourceb)
 if ismember(lower(sourceb(i)),alphabet)
 sigb=[sigb sourceb(i)];
 end
 end
 sigb;
 pithanothta=my_posibilities(sigb);
 alpha=my_epektash(alphabet);
 pr=my_prob(pithanothta);
 dict=mey_huffmandict(alpha,pr);
 sigbb={};
 for i=1 : length(sourceb)-1
 if ismember(lower(sourceb(i)),alphabet) &&
ismember(lower(sourceb(i+1)),alphabet)
 sigbb=[sigbb strcat(sourceb(i),sourceb(i+1))];
 end
 end
 encob=myhuffmanenco(sigbb,dict);
 dsigb = my_huffmandeco(encob,dict);
 entropia=0;
 for l=1 : 676
 entropia= entropia +pr(l)*log2(1/pr(l));
 end
 mesos={};

 for j=1 :676
 mesos(j) = dict(j,2);

 end
 length=cellfun('length',mesos);
 true_mesos=0;
 for k=1 : 676
 true_mesos= true_mesos + pr(1,k)*length(1,k)';
 end

 apodosh=entropia/true_mesos;
function [dict,avglen] = mey_huffmandict(sym, prob)
% Check an einai cell array an oxi to metatrepw
if ~iscell(sym)
 sym = num2cell(sym);
end
% Check an to prob einai opws prepei
validateattributes(prob, {'double'}, ...
 {'real','vector','nonnegative','<=',1}, ...
 'my_huffmandict', 'prob');
% Arxikopoihsh tou n_ary sto 2 afou tha asxolithoume mono me 2adiki
% kodikopoihsh
n_ary = 2;
% Efoson asxoloumaste me 2adiki kodikopoihsh mono to variance ginete
%set se
% max
 variance = 'max';
%Dhmiourgoume dedro me ta symbols kai tis adistixes pithanotites
h_tree = struct('signal', [], 'probability', [],...
 'child', [], 'code', [], 'origOrder', -1);
for i = 1:length( sym )
 h_tree(i).signal = sym{i};
 h_tree(i).probability = prob(i); 
 h_tree(i).origOrder = i;
end
% Sort ta symbols kai tis pithanotites basismenoi stis pithanothtes
[~, i] = sort(prob);
h_tree = h_tree(i);
% Dhmiourgia dedrou
h_tree = create_tree(h_tree, n_ary, variance);
[~,dict,avglen] = create_dict(h_tree,{},0);
% Kanoume sort to dict
[~, dictsortorder] = sort([dict{:,4}]);
lenDict = length(dictsortorder);
finaldict = cell(lenDict, 2);
for i=1:length(dictsortorder)
 finaldict{i,1} = dict{dictsortorder(i), 1};
 finaldict{i,2} = dict{dictsortorder(i), 2};
end
dict = finaldict;
%Dhmiourgoume dedro me metablites n_ary=2 kai variance=max
function h_tree = create_tree(h_tree, n_ary, variance)
% Testaroume an to mhkos tou dedrou einai 1 auto shmenei oti iparxei
%1
% kombos sto array twn kombwn kai to programma termatizei
numRemNodes = length(h_tree);
if( numRemNodes <= 1)
 return;
end
%Sindiazontas tous 2 kombous me tis mikroteres pithanotites se
enan,enw
%tautoxrona aferoume tous 2 sindazomenous kombous kai prosthatoume
ton neo.
%Episis anathetoume diadikes times gia na prokipsei i kodikopoihsh
temp = struct('signal', [], 'probability', 0, ...
 'child', [], 'code', []);
numNodesToComb = rem(numRemNodes-1, n_ary-1) + 1;
if numNodesToComb == 1
 numNodesToComb = n_ary;
end
for i = 1:numNodesToComb
 if isempty(h_tree), break; end
 temp.probability = temp.probability + h_tree(1).probability;
 temp.child{i} = h_tree(1);
 temp.origOrder = -1;
 h_tree(1) = [];
end
if( strcmpi(variance, 'min') == 1 )
 h_tree = insertMinVar(h_tree, temp);
else
 h_tree = insertMaxVar(h_tree, temp);
end
%Dhmiourgoume ena neo dedro me tous idh meiwmenous komvous
h_tree = create_tree(h_tree, n_ary, variance);
end
%Auth h synarthsh eisagei ton kombo sthn sortarismenh lista. An
uparxei idh
%kombos me tin idia pithanotita emfanisis me ton idio kombo tote o
%kenurgious kombos topothetite meta ton palio.
function h_tree = insertMaxVar(h_tree, newNode)
i = 1;
while i <= length(h_tree) && ...
 newNode.probability > h_tree(i).probability
 i = i+1;
end
h_tree = [h_tree(1:i-1) newNode h_tree(i:end)];
end
%Se auth th synartisi ginete oti ginete stin apopanw me thn diafora
oti o
%neos kombos benei prin apo ton palio.
function h_tree = insertMinVar(h_tree, newNode)
i = 1;
while i <= length(h_tree) && ...
 newNode.probability >= h_tree(i).probability
 i = i+1;
end
h_tree = [h_tree(1:i-1) newNode h_tree(i:end)];
end
%Edw dimiourgite o kodikas huffman gia tous kombous tou dedrou
ksekinodas
%apo ta fila tou dedrou
function [h_tree,dict,b] = create_dict(h_tree,dict,b)
%Tsekaroume an o kombos einai fulo tou dedrou an einai tote tou
bazoume ta
%adistixa 0,1
if isempty(h_tree.child)
 dict{end+1,1} = h_tree.signal;
 dict{end, 2} = h_tree.code;
 dict{end, 3} = length(h_tree.code);
 dict{end, 4} = h_tree.origOrder;
 b = b + length(h_tree.code)*h_tree.probability;
 return;
end
num_childrens = length(h_tree.child);
for i = 1:num_childrens
 h_tree.child{i}.code = [h_tree(end).code, (num_childrens-i)];
 [h_tree.child{i}, dict, b] = ...
 create_dict(h_tree.child{i}, dict, b);
end
end
end
function enco = myhuffmanenco(text, dict)
enco = [];
for i=1:length(text)

 cc = text(i);
 % Ftixanoume to t me deiktes.

 t = strcmp(cc, dict);
 for j=1:length(t)
 if (t(j,1) == 1)
 t(j,1) = 0;
 t(j,2) = 1;
 end
 end
 % Briskoume ton katallhlo deikti
 enco = [enco dict{t}];
end
end
function text = my_huffmandeco(enco, dict)
i = 1;
enco = sprintf('%d', enco);
while i <= length(dict)
 % Metatroph eisodou se string
 code = sprintf('%d', dict{i,2});

 patt = ['^' code];
 if regexp(enco, patt)
 try
 text = [text {dict{i,1}}];
 catch
 text = {dict{i,1}};
 end
 % Arxikopoihsh
 i = 1;

 le = length(enco);
 lc = length(code) + 1;
 text=cell2mat(text);

 if lc > le
 break;
 else
 enco = enco(lc:le);
 end
 else
 i = i + 1;
 end

end
text=text';
end
function [posib] = my_posibilities(sigb)
sigmat=cell2mat(sigb);
counta=0;
countb=0;
countc=0;
countd=0;
counte=0;
countf=0;
countg=0;
counth=0;
counti=0;
countj=0;
countk=0;
countl=0;
countm=0;
countn=0;
counto=0;
countp=0;
countq=0;
countr=0;
counts=0;
countt=0;
countu=0;
countv=0;
countw=0;
countx=0;
county=0;
countz=0;
for i=1 : length(sigmat)
if sigmat(i) == 'a'
 counta=counta+1;
elseif sigmat(i) == 'b'
 countb=countb+1;
elseif sigmat(i) == 'c'
 countc=countc+1;
elseif sigmat(i) == 'd'
 countd=countd+1;
elseif sigmat(i) == 'e'
 counte=counte+1;

elseif sigmat(i) == 'f'
 countf=countf+1;

elseif sigmat(i) == 'g'
 countg=countg+1;
elseif sigmat(i) == 'h'
 counth=counth+1;
elseif sigmat(i) == 'i'
 counti=counti+1;

elseif sigmat(i) == 'j'
 countj=countj+1;

elseif sigmat(i) == 'k'
 countk=countk+1;

elseif sigmat(i) == 'l'
 countl=countl+1;
 
elseif sigmat(i) == 'm'
 countm=countm+1;

elseif sigmat(i) == 'n'
 countn=countn+1;

elseif sigmat(i) == 'o'
 counto=counto+1;

elseif sigmat(i) == 'p'
 countp=countp+1;

elseif sigmat(i) == 'q'
 countq=countq+1;

elseif sigmat(i) == 'r'
 countr=countr+1;

elseif sigmat(i) == 's'
 counts=counts+1;

elseif sigmat(i) == 't'
 countt=countt+1;

elseif sigmat(i) == 'u'
 countu=countu+1;

elseif sigmat(i) == 'v'
 countv=countv+1;

elseif sigmat(i) == 'w'
 countw=countw+1;
elseif sigmat(i) == 'x'
 countx=countx+1;

elseif sigmat(i) == 'y'
 county=county+1;

elseif sigmat(i) == 'z'
 countz=countz+1;

end
end
posib(1)=counta/length(sigb);
posib(2)=countb/length(sigb);
posib(3)=countc/length(sigb);
posib(4)=countd/length(sigb);
posib(5)=counte/length(sigb);
posib(6)=countf/length(sigb);
posib(7)=countg/length(sigb);
posib(8)=counth/length(sigb);
posib(9)=counti/length(sigb);
posib(10)=countj/length(sigb);
posib(11)=countk/length(sigb);
posib(12)=countl/length(sigb);
posib(13)=countm/length(sigb);
posib(14)=countn/length(sigb);
posib(15)=counto/length(sigb);
posib(16)=countp/length(sigb);
posib(17)=countq/length(sigb);
posib(18)=countr/length(sigb);
posib(19)=counts/length(sigb);
posib(20)=countt/length(sigb);
posib(21)=countu/length(sigb);
posib(22)=countv/length(sigb);
posib(23)=countw/length(sigb);
posib(24)=countx/length(sigb);
posib(25)=county/length(sigb);
posib(26)=countz/length(sigb);
end
function [alpha] = my_epektash(alphabet)
d={};
for i=1 : length(alphabet)
for j=1 : 26
 d = [d strcat(alphabet(i),alphabet(j))];
end
end
alpha=d;
end
function [prob] = my_prob(p)
d=[];
for i=1 : length(p)
for j=1 : 26
 d = [d p(i)*p(j)];
end
end
prob=d;
end