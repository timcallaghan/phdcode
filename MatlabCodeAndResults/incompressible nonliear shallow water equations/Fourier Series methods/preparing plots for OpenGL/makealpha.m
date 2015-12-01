function Achannel = makealpha(PNG,toll)
%MAKEALPHA Given a PNG file with no alpha channel
% this script constucts an alpha channel by examining
% the RGB values of each pixel...if it's white it makes it
% transparent...else it's made opaque

picsize=length(PNG);
Achannel=zeros(picsize);

for i=1:picsize
    for j=1:picsize
        if (PNG(i,j,1)>=toll & PNG(i,j,2)>=toll & PNG(i,j,3)>=toll)
            Achannel(i,j)=1.0;
        else
            Achannel(i,j)=0.0;
        end
    end
end
