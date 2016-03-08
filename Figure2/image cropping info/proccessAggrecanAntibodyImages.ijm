

open("E:\\2015 Data\\2015-01-29_MSCAggAb --Don8and9\\2015-01-29_donor9_day7CMplus\\examples\\14\\gfp014.tif");
run("Duplicate...", "title=gfp014-16.tif duplicate range=16-16");
selectWindow("gfp014.tif");
close();

open("E:\\2015 Data\\2015-01-29_MSCAggAb --Don8and9\\2015-01-29_donor9_day7CMplus\\examples\\14\\dapi014.tif");
run("Duplicate...", "title=dapi014-16.tif duplicate range=16-16");
selectWindow("dapi014.tif");
close();


saveAs("Text", "C:\\Program Files (x86)\\ImageJ\\macros\\proccessAggrecanAntibodyImages.ijm");
selectWindow("dapi014-16.tif");
selectWindow("gfp014-16.tif");


run("Images to Stack", "name=Stack title=[] use keep");
run("Make Composite", "display=Composite");
selectWindow("Stack");
run("Grays");
run("Magenta");
//run("Brightness/Contrast...");

setMinAndMax(810, 2063);
call("ij.ImagePlus.setDefault16bitRange", 0);

setMinAndMax(2271, 22987);
call("ij.ImagePlus.setDefault16bitRange", 0);

//setTool("rectangle");
makeRectangle(196, 92, 700, 700);
run("Crop");
makeRectangle(58, 4, 472, 472);
run("Crop");



//setTool("rectangle");
makeRectangle(196, 92, 700, 700);
run("Crop");
makeRectangle(58, 4, 472, 472);
run("Crop");


run("Duplicate...", "title=trans014-1.tif duplicate range=7-7");
saveAs("PNG", "E:\\2015 Data\\2015-01-29_MSCAggAb --Don8and9\\2015-01-29_donor9_day7CMplus\\examples\\14\\trans014-1_cropped.png");
selectWindow("trans014.tif");


