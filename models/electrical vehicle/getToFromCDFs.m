function [ToCDF,FromCDF] = getToFromCDFs(ToPDF,FromPDF)

ToCDF = zeros(1,length(ToPDF));
FromCDF = zeros(1,length(FromPDF));
ToNormalize = sum(ToPDF);
FromNormalize = sum(FromPDF);
ToPDF = ToPDF/ToNormalize;
FromPDF = FromPDF/FromNormalize;

ToCDF(1) = ToPDF(1);
FromCDF(1) = FromPDF(1);
for i=2:1:length(ToPDF)
    ToCDF(i) = ToCDF(i-1)+ToPDF(i);
    FromCDF(i) = FromCDF(i-1)+FromPDF(i);
end