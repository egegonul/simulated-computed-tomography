
function [filtered_image,unfiltered_image]=back_project(load_type,projections_file,image_size)

%depending on the data loading type load the data
%load_type==0 for .mat loading
%load_type==1 for .txt loading
    if ~load_type
        load(projections_file,'projections');
        projections=projections';
        num_beams=size(projections,2);
        num_proj=size(projections,1);
        projections=projections';
    else
        fileID = fopen(projections_file,'r');
         formatSpec = '%f';
         A = fscanf(fileID,formatSpec);
        num_proj=A(1);
        num_beams=A(2);
        projections=reshape(A(3:end),[num_beams+1,num_proj]);
        projections=projections(2:end,:);
    end

    %initialize the backprojected image as zeros
    filtered_image=zeros(image_size,image_size);
    unfiltered_image=zeros(image_size,image_size);

    
    step_size=180/num_proj;
    M=size(filtered_image,1);
    beams=linspace(-M/sqrt(2),M/sqrt(2),num_beams);
    thetas=1:step_size:180;
    all_data=cell(length(thetas),1);


    %recalculate adress and distance values in case projections come from 
    %outside
    for ww=1:length(thetas)
        theta=thetas(ww);
        data=cell(num_beams,3);
        k=0;
        for t=beams
            k=k+1;
            %define x and y vectors
            x=(-M/2:1:M/2)';
            y=(-M/2:1:M/2)';
        
            %find intersections points for x and y
            yx=(t-x*cosd(theta))/(sind(theta)+0.00001);
            inds=yx<M/2&yx>-M/2;
            yx(abs(yx)<0.001)=0;
            yx=inds.*yx;
            inds=find(yx);
            relevant_x=[x(inds) yx(inds)];
            
            xy=(t-y*sind(theta))/(cosd(theta)+0.00001);
            xy(abs(xy)<0.001)=0;
            inds=xy<M/2&xy>-M/2;
            xy=inds.*xy;
            inds=find(xy);
            relevant_y=[xy(inds) y(inds)];
        
            %reduce the intersection points to relevant points
            relevant_points=unique([relevant_x;relevant_y],'rows');
            relevant_points=sortrows(relevant_points);
          
            %calculate midpoints and distances only if there are nonzero relevant
            %points
            if(size(relevant_points,1)>1)
                %calculate distance and midpoint
                dist=[];
                mid=[];
                for i=1:size(relevant_points,1)-1
                    x_cur=relevant_points(i,1);
                    y_cur=relevant_points(i,2);
                    x_next=relevant_points(i+1,1);
                    y_next=relevant_points(i+1,2);
                    dist=[dist; sqrt((x_next-x_cur)^2+(y_next-y_cur)^2)];
                    mid=[mid;(x_cur+x_next)/2,(y_cur+y_next)/2];
                end
                data{k,2}=dist;
                
                
                %detect the adresses
                adress=[M/2-floor(mid(:,2)), M/2+ceil(mid(:,1))];
                data{k,3}=adress;
       
            end
    
        all_data{ww}=data;
    end
    end


    %filter all the projections
    filtered_proj=ram_lak_filter(projections,num_proj);

    % back_project without filtering
    for ww=1:(length(all_data))
        data=all_data{ww};
        for i=1:num_beams
            proj=projections(i,ww);
            if ~proj|~isempty(proj)
                dist=data{i,2};
                adress=data{i,3};
                for j=1:length(dist)
                    unfiltered_image(adress(j,1),adress(j,2))=unfiltered_image(adress(j,1),adress(j,2))+dist(j)*proj;
                end
            end
        end
    end
    

    %back_project with filtering
    for ww=1:(length(all_data))
        data=all_data{ww};
        for i=1:num_beams
            proj=filtered_proj(i,ww);
            if ~proj|~isempty(proj)
                dist=data{i,2};
                adress=data{i,3};
                for j=1:length(dist)
                    filtered_image(adress(j,1),adress(j,2))=filtered_image(adress(j,1),adress(j,2))+dist(j)*proj;
                end
            end
        end
    end
    
  
    subplot(2,1,1);
    nrm=max(filtered_image,[],'all');
    filtered_image=filtered_image./nrm;
    imshow(filtered_image);
    title("filtered backprojection")
    subplot(2,1,2);
    nrm=max(unfiltered_image,[],'all');
    unfiltered_image=unfiltered_image./nrm;
    imshow(unfiltered_image);
    title("unfiltered backprojection")

    
end
    



