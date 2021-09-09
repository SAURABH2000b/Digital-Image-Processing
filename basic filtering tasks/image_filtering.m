box_filter_width=5;
tent_filter_width=5;
gaussian_filter_width=5;
B_Spline_cubic_filter_width=5;
catmull_rom_cubic_filter_width=5;
mitchell_netravali_cubic_filter_width=5;

sampling_rate=3;

box_filter = ones(box_filter_width, box_filter_width)*(1/9);

tent_equation = @(x,y,width) (-4*abs(x)/(width*width)+2/width)*(-4*abs(y)/(width*width)+2/width);
tent_filter = ones(tent_filter_width, tent_filter_width);
for i=1:tent_filter_width
    for j=1:tent_filter_width
        xx=j-ceil(tent_filter_width/2);
        yy=i-ceil(tent_filter_width/2);
        tent_filter(i,j)=tent_equation(xx,yy,tent_filter_width);
    end
end

gaussian_equation = @(x,y) (1/(2*pi)*exp(-1*(x*x+y*y)/2));
gaussian_filter = ones(gaussian_filter_width, gaussian_filter_width);
for i=1:gaussian_filter_width
    for j=1:gaussian_filter_width
        xx=(j-ceil(gaussian_filter_width/2))/sampling_rate;
        yy=(i-ceil(gaussian_filter_width/2))/sampling_rate;
        gaussian_filter(i,j)=gaussian_equation(xx,yy);
    end
end

B_Spline_cubic_equation_1 = @(x) ((1/6)*(-3*(1-abs(x))*(1-abs(x))*(1-abs(x)) + 3*(1-abs(x))*(1-abs(x)) + 3*(1-abs(x)) + 1));
B_Spline_cubic_equation_2 = @(x) ((1/6)*(2-abs(x))*(2-abs(x))*(2-abs(x)));
B_Spline_filter = ones(B_Spline_cubic_filter_width, B_Spline_cubic_filter_width);
for i=1:B_Spline_cubic_filter_width
    for j=1:B_Spline_cubic_filter_width
        xx=(j-ceil(B_Spline_cubic_filter_width/2))/sampling_rate;
        yy=(i-ceil(B_Spline_cubic_filter_width/2))/sampling_rate;
        if abs(xx)>2 || abs(yy)>2
            B_Spline_filter(i,j)=0;
            continue;
        end
        if xx>=-1 && xx<=1
            if yy>=-1 && yy<=1
                B_Spline_filter(i,j)=B_Spline_cubic_equation_1(xx)*B_Spline_cubic_equation_1(yy);
            else
                B_Spline_filter(i,j)=B_Spline_cubic_equation_1(xx)*B_Spline_cubic_equation_2(yy);
            end
        else
            if yy>=-1 && yy<=1
                B_Spline_filter(i,j)=B_Spline_cubic_equation_2(xx)*B_Spline_cubic_equation_1(yy);
            else
                B_Spline_filter(i,j)=B_Spline_cubic_equation_2(xx)*B_Spline_cubic_equation_2(yy);
            end
        end
    end
end

Catmull_Rom_cubic_equation_1 = @(x) ((1/2)*(-3*(1-abs(x))*(1-abs(x))*(1-abs(x)) + 4*(1-abs(x))*(1-abs(x)) + 1-abs(x)));
Catmull_Rom_cubic_equation_2 = @(x) ((1/2)*((2-abs(x))*(2-abs(x))*(2-abs(x))-(2-abs(x))*(2-abs(x))));
Catmull_Rom_cubic_filter= ones(catmull_rom_cubic_filter_width, catmull_rom_cubic_filter_width);
for i=1:catmull_rom_cubic_filter_width
    for j=1:catmull_rom_cubic_filter_width
        xx=(j-ceil(catmull_rom_cubic_filter_width/2))/sampling_rate;
        yy=(i-ceil(catmull_rom_cubic_filter_width/2))/sampling_rate;
        if abs(xx)>2 || abs(yy)>2
            Catmull_Rom_cubic_filter(i,j)=0;
            continue;
        end
        if xx>=-1 && xx<=1
            if yy>=-1 && yy<=1
                Catmull_Rom_cubic_filter(i,j)=Catmull_Rom_cubic_equation_1(xx)*Catmull_Rom_cubic_equation_1(yy);
            else
                Catmull_Rom_cubic_filter(i,j)=Catmull_Rom_cubic_equation_1(xx)*Catmull_Rom_cubic_equation_2(yy);
            end
        else
            if yy>=-1 && yy<=1
                Catmull_Rom_cubic_filter(i,j)=Catmull_Rom_cubic_equation_2(xx)*B_Spline_cubic_equation_1(yy);
            else
                Catmull_Rom_cubic_filter(i,j)=Catmull_Rom_cubic_equation_2(xx)*Catmull_Rom_cubic_equation_2(yy);
            end
        end
    end
end

Mitchell_Netravali_cubic_filter_equation_1 = @(x) ((1/18)*(-21*(1-abs(x))*(1-abs(x))*(1-abs(x)) + 27*(1-abs(x))*(1-abs(x)) + 9*(1-abs(x)) + 1));
Mitchell_Netravali_cubic_filter_equation_2 = @(x) ((1/18)*(7*(2-abs(x))*(2-abs(x))*(2-abs(x))-6*(2-abs(x))*(2-abs(x))));
Mitchell_Netravali_cubic_filter = ones(mitchell_netravali_cubic_filter_width, mitchell_netravali_cubic_filter_width);
for i=1:mitchell_netravali_cubic_filter_width
    for j=1:mitchell_netravali_cubic_filter_width
        xx=(j-ceil(mitchell_netravali_cubic_filter_width/2))/sampling_rate;
        yy=(i-ceil(mitchell_netravali_cubic_filter_width/2))/sampling_rate;
        if abs(xx)>2 || abs(yy)>2
            Mitchell_Netravali_cubic_filter(i,j)=0;
            continue;
        end
        if xx>=-1 && xx<=1
            if yy>=-1 && yy<=1
                Mitchell_Netravali_cubic_filter(i,j)=Mitchell_Netravali_cubic_filter_equation_1(xx)*Mitchell_Netravali_cubic_filter_equation_1(yy);
            else
                Mitchell_Netravali_cubic_filter(i,j)=Mitchell_Netravali_cubic_filter_equation_1(xx)*Mitchell_Netravali_cubic_filter_equation_2(yy);
            end
        else
            if yy>=-1 && yy<=1
                Mitchell_Netravali_cubic_filter(i,j)=Mitchell_Netravali_cubic_filter_equation_2(xx)*Mitchell_Netravali_cubic_filter_equation_1(yy);
            else
                Mitchell_Netravali_cubic_filter(i,j)=Mitchell_Netravali_cubic_filter_equation_2(xx)*Mitchell_Netravali_cubic_filter_equation_2(yy);
            end
        end
    end
end

box_filter = box_filter/sum(box_filter,'all');
tent_filter = tent_filter/sum(tent_filter,'all');
gaussian_filter = gaussian_filter/sum(gaussian_filter,'all');
B_Spline_filter = B_Spline_filter/sum(B_Spline_filter,'all');
Catmull_Rom_cubic_filter = Catmull_Rom_cubic_filter/sum(Catmull_Rom_cubic_filter,'all');
Mitchell_Netravali_cubic_filter = Mitchell_Netravali_cubic_filter/sum(Mitchell_Netravali_cubic_filter,'all');

test_image = imread("sample_image_2.jfif");
test_image = double(test_image);
padding = floor(box_filter_width/2);
padded_image = pad(test_image,padding);
filtered_image_BF = filtering(padded_image, box_filter, padding, box_filter_width);
filtered_image_TF = filtering(padded_image, tent_filter, padding, tent_filter_width);
filtered_image_GF = filtering(padded_image, gaussian_filter, padding, gaussian_filter_width);
filtered_image_BSF = filtering(padded_image, B_Spline_filter, padding, B_Spline_cubic_filter_width);
filtered_image_CRF = filtering(padded_image, Catmull_Rom_cubic_filter, padding, catmull_rom_cubic_filter_width);
filtered_image_MNF = filtering(padded_image, Mitchell_Netravali_cubic_filter, padding, mitchell_netravali_cubic_filter_width);

figure, imshow(filtered_image_BF);
figure, imshow(filtered_image_TF);
figure, imshow(filtered_image_GF);
figure, imshow(filtered_image_BSF);
figure, imshow(filtered_image_CRF);
figure, imshow(filtered_image_MNF);

function pad_image = pad(img, padding)
image_dimensions = size(img);
pad_image_dimx = image_dimensions(1) + 2*padding;
pad_image_dimy = image_dimensions(2) + 2*padding;
pad_image = zeros(pad_image_dimx, pad_image_dimy,3);
pad_image(padding+1:end-padding, padding+1:end-padding,:) = img;
end

function filtered_image = filtering(img, filter, padding, filter_width)
image_dimensions = size(img);
image_dimx = image_dimensions(1);
image_dimy = image_dimensions(2);
filtered_image = zeros(image_dimx-2*padding, image_dimy-2*padding);
for i=1:(image_dimx-2*padding)
    for j=1:(image_dimy-2*padding)
        for k=1:3
            filtered_image(i, j, k) = sum(filter.*img(i:i+filter_width-1,j:j+filter_width-1, k), 'all');
        end
    end
end
filtered_image = uint8(filtered_image);
end
