function cout = rename_clustering(cin)


current = 1;
n = numel(cin);
num = max(cin);

map = zeros(num,1);
map(1) = cin(1);
seen = [cin(1)];
next = 2;
for i = 2:n
    if ismember(cin(i),seen)
        continue
    else
        map(next) = cin(i);
        next = next+1;
        seen = [seen;cin(i)];
    end
end

cout = zeros(n,1);
for j = 1:num
    set = find(cin == map(j));
    cout(set) = j;
end

end