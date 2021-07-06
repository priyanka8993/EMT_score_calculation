function [A5, A6] = MLR3(A1, A2, A3, A4)

[a b]=size(A4);
if length(A4(:,1))==length(A2(:,1))
    A5=mnrval(A1(:,1),A4','model','ordinal');
    A6=zeros(b,1);
    for j=1:b
        if isequaln(A5(j,:),NaN.*ones(1,length(A1(:,1))-1))==1
            A6(j,1)=NaN;
        else
            A6(j,1)=find(A5(j,:)==max(A5(j,:)));
        end
    end
else
    'Hello'
    for j=1:length(A2(:,1))
        temp=find(strcmp(A2(j,1),A3)==1);
        index(j)=temp(1);
    end
    A7=zeros(length(A2(:,1)),b);
    for j=1:length(A2(:,1))
        A7(j,:)=A4(index(j),:);
    end
    A5=mnrval(A1(:,1),A7','model','ordinal');
    A6=zeros(b,1);
    for j=1:b
        if isequaln(A5(j,:),NaN.*ones(1,length(A1(:,1))-1))==1
            A6(j,1)=NaN;
        else
            A6(j,1)=find(A5(j,:)==max(A5(j,:)));
        end
    end
end

end