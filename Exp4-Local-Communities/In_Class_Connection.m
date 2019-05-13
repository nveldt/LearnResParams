%% Go back to the basics and compute the volumes this time, not just the set sizes
fid = fopen('Facebook_Sets.txt','r');
addpath('~/data/Facebook100/')
f = fgets(fid);
names = split(f);
minSize = 20;
YearSets = zeros(700,7);
next = 1;


good = 0;
for i = 1:numel(names)
    
d2005 = [];
d2006 = [];
d2007 = [];
d2008 = [];
d2009 = [];
din2005 = [];
din2006 = [];
din2007 = [];
din2008 = [];
din2009 = [];
    name = names{i};                            
    load(strcat(name))
    volA = sum(nonzeros(A));
    d = full(sum(A,1)');
    
    
    yearlabels = 2005:2009;
    numyears = numel(yearlabels);
    
    
    for j = 1:numel(yearlabels)
        
        y = yearlabels(j);
        year = find(local_info(:,6) == y);
        Ayear = A(year,year);
        
        % Get the in-year degree
        dyear = full(sum(Ayear,1)');
        
        if y == 2005
            d2005 = [d2005; d(year)];
        elseif y == 2006
            d2006 = [d2006; d(year)];
        elseif y == 2007
            d2007 = [d2007; d(year)];
        elseif y == 2008
            d2008 = [d2008; d(year)];
        else
            d2009 = [d2009; d(year)];
        end

        if y == 2005
            din2005 = [din2005; dyear];
        elseif y == 2006
            din2006 = [din2006; dyear];
        elseif y == 2007
            din2007 = [din2007; dyear];
        elseif y == 2008
            din2008 = [din2008; dyear];
        else
            din2009 = [din2009; dyear];
        end
    end

rat2006 = mean(din2006)/mean(d2006);
rat2007 = mean(din2007)/mean(d2007);
rat2008 = mean(din2008)/mean(d2008);
rat2009 = mean(din2009)/mean(d2009);

if rat2006 < rat2007 && rat2007 < rat2008 && rat2008 < rat2009
    good = good+1;
    
    fprintf('\n%s\n',name)
    fprintf('2006: %f, %f, %f \n',mean(d2006),mean(din2006),mean(din2006)/mean(d2006))
    fprintf('2007: %f, %f, %f \n',mean(d2007),mean(din2007),mean(din2007)/mean(d2007))
    fprintf('2008: %f, %f, %f \n',mean(d2008),mean(din2008),mean(din2008)/mean(d2008))
    fprintf('2009: %f, %f, %f \n',mean(d2009),mean(din2009),mean(din2009)/mean(d2009))

else
    % Special print for the exceptions where the average within-class
    % friendship ratio does not follow a monotonic pattern across
    % graduation years
    fprintf('\n%s: NONSTANDARD\n',name)
    fprintf('2006: %f, %f, %f \n',mean(d2006),mean(din2006),mean(din2006)/mean(d2006))
    fprintf('2007: %f, %f, %f \n',mean(d2007),mean(din2007),mean(din2007)/mean(d2007))
    fprintf('2008: %f, %f, %f \n',mean(d2008),mean(din2008),mean(din2008)/mean(d2008))
    fprintf('2009: %f, %f, %f \n',mean(d2009),mean(din2009),mean(din2009)/mean(d2009))
end

end

