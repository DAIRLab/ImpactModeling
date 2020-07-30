%%make updates to data set

load('squareDataPhilUpdated.mat');
threshold = 4;
count = 0;
for i = 1:length(squareDataPhilUpdated)
    disp(squareDataPhilUpdated(i).states(12))
    if abs(squareDataPhilUpdated(i).states(12)) > threshold && abs(squareDataPhilUpdated(i).states(6)) < 30 
        count = count + 1;
        squareDataPhilChanged(count) = squareDataPhilUpdated(i);
    end
end

