Np=5; %Number of interp points along one direction in one elem
K=.1*40; %Non-zero elements in one direction
N=K^2; %Total non-zero elements

%Stencil for particular point with self ref and num flux ref
%To right
Choose(1,:)=1:Np*2;
%To left
Choose(2,:)=-(Np-1):Np;
%To right wrap
Choose(3,:)=[1:Np,(1:Np)-((K-1)*Np)];
%To left wrap
Choose(4,:)=[(1:Np)+((K-1)*Np),1:Np];

%Random matrix vals, would be stiffness-num flux etc
s=rand(Np^2*2*Np*N,1);
%Location of nz's in sparse global stencil
%2*Np refs per point on strip (self ref with forward or backward ref for num flux)
%Np^2 points per elem, with N elems
i=reshape(repmat(1:Np^2*N,2*Np,1),[],1);

%Ditto from 'i'
j=zeros(Np^2*N,Np*2);

%Somewhat involved provess, blech
%Iterate by elements in the y, for each clump of elems at a particular y,
%go down each "stream" of strips at a particular Np_y in that element clump
%Within each element_x's portion of that stream, step through each interp
%point. If the element_x of a stream is on the left or right, it may wrap to
%the right or left respectively for periodic BCs on that stream
for elem_y=(0:K-1)*Np*K*Np
    for stream_y=(0:Np-1)*Np*K
        elem_x=0;
        for point=1:Np
            j(elem_y+stream_y+elem_x+point,:)=Choose(randi([0,1],1)*3+1,:);
        end
        for elem_x=(1:K-2)*Np
            for point=1:Np
                j(elem_y+stream_y+elem_x+point,:)=Choose(randi(2,1),:);
            end
        end
        elem_x=(K-1)*Np;
        for point=1:Np
            j(elem_y+stream_y+elem_x+point,:)=Choose(randi([2,3],1),:);
        end
    end
end

%Enforce first Np points in top left strip don't wrap to avoid messing with
%the global matrix
j(1:Np,:)=repmat(Choose(1,:),Np,1);
%Ditto for bottom right strip
j(end-(Np-1):end,:)=repmat(Choose(2,:),Np,1);

advance=reshape(repmat(0:Np:K*(Np*(Np*K))-1,Np,1),[],1);
j=reshape(bsxfun(@plus,j,advance)',[],1);

A=sparse(i,j,s,Np^2*N,Np^2*N);


