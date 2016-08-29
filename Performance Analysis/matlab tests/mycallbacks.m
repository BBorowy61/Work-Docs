function a = mycallbacks(str)
a =str2func(str);

function myclick(varargin)
disp('Single click function')


function my2click(varargin)
disp('Double click function')


function mymoused(varargin)
disp('Tou have reached the mouse down finction')
disp('The X position is: ')
double(varargin{5})
disp('The Y position is: ')
double(varargin{6})



