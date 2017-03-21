function [Away,AwayStart] = getAwayPeriod(ToCDF,FromCDF)

Away = zeros(1,length(ToCDF));
% Pick from distribution
eventTime1 = find(rand<ToCDF,1,'first');
eventTime2 = find(rand<FromCDF,1,'first');
AwayStart = min(eventTime1,eventTime2);
% Away hours
for i = AwayStart:1:max(eventTime1,eventTime2)
    Away(i) = true;
end