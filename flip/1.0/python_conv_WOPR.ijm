dir = getDirectory("Choose a Directory: ");

// used for local testing
// dir = "C:/Users/David.Moller/Desktop/2019-07-12/";
// dir = "C:/Users/David.Moller/Desktop/2019-07-12/2019-07-12__20-38-24-656";

//print("start")

setBatchMode(true);
count = 0;
countFiles(dir);
n = 0;
processFiles(dir);

//print(count + " files processed");


//****************************************************************************************************
function countFiles(dir) {

    list = getFileList(dir);

    for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/")) {
           // print(list[i]);
            countFiles(""+dir+list[i]);
        }
        else
            count++;
    }
    //print(count);
}

function processFiles(dir) {
    list1 = getFileList(dir);

    for (i=0; i<list1.length; i++) {
        if (endsWith(list1[i], "/")) {
        //print(list[i]);
        processFiles(""+dir+list1[i]);
        }
        else {
            showProgress(n++, count);
            path = dir+list1[i];
            processFile(path);
        }
    //print("Looking for Bin_" + list1[1])
    }
}

// code not needed as the python is handling the conversion
//function processFile(path) {
//
//    print ("Looking for BIN file")
//
//    file=split(path, ".");
//    if (endsWith(path, ".bin")) {
//        run("Raw...", 'open=[path] image=[8-bit] width=1936 height=1216 offset=0 number=1 gap=0');
//        saveAs("png", file[0] + ".PNG");
//        close();
//        print(file[0])
//        print("found Bin file")
//    }
//}

//*****************************************************************************************************************

function PS2(dir) {
    //print("in the PS2")

    list = getFileList(dir);

    for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/"))
            countFiles(""+dir+list[i]);
        else
            count++;
    }
}

function processFiles(dir) {
    //print("in the processFiles again?");
    list = getFileList(dir);

    for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/"))
            processFiles(""+dir+list[i]);
        else {
            showProgress(n++, count);
            path = dir+list[i];
            processFile(path);
        }
    }
    
    file  = dir;
    file2 = replace(file, "\\", "@");
    file3 = split(file2, "@");
    //print(file2);
    //print(file3.length);
    if (file3.length == 7) {
        file4 = split(file3[6], "/");
		print(file4[0]);
        selectWindow("Results");
        saveAs("Results", dir + file4[0] + ".csv");
        close("Results");
        call("java.lang.System.gc");
        //print("folder processed properly");
    }
    else {
        //print("folder skipped because array out of bounds");
    }
}

function processFile(path) {

    if (endsWith(path, ".png")) {
        open(path);
        makeRectangle(0,405,1936,405);
        run("Crop");
        setThreshold(0, 7);
        run("Create Selection");
        run("Measure");
        run("Select None");
        setThreshold(8, 10);
        run("Create Selection");
        run("Measure");
        run("Select None");
        setThreshold(11, 14);
        run("Create Selection");
        run("Measure");
        run("Select None");
        setThreshold(15, 19);
        run("Create Selection");
        run("Measure");
        run("Select None");
        setThreshold(20, 255);
        run("Create Selection");
        run("Measure");
        run("Select None");
        close();	
    }
}



