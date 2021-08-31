package com.example.vulkantriangle;

import androidx.appcompat.app.AppCompatActivity;
import android.app.NativeActivity;

import android.content.Intent;
import android.os.Bundle;
import android.widget.Toast;

import static android.Manifest.permission.READ_EXTERNAL_STORAGE;
import static android.Manifest.permission.WRITE_EXTERNAL_STORAGE;
//import com.example.vulkanTriangle.R;

import java.io.File;

public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        System.out.println("Working Directory = " +
                System.getProperty("user.dir"));
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

//        if (loadNativeLibrary(getResources().getString(R.string.native_lib_name))) {
//            File externalFilesDir = getExternalFilesDir("");
//            if (externalFilesDir != null) {
//                initFilePath(externalFilesDir.toString());
//            }
//        }

        requestPermissions(new String[]{ WRITE_EXTERNAL_STORAGE, READ_EXTERNAL_STORAGE}, 1);

//        Intent intent = new Intent(MainActivity.this, NativeActivity.class);
//        startActivity(intent);
//        finishAffinity();
    }

//    private boolean loadNativeLibrary(String nativeLibName) {
//        boolean status = true;
//
//        try {
//            System.loadLibrary(nativeLibName);
//        } catch (UnsatisfiedLinkError e) {
//            Toast.makeText(getApplicationContext(), "Native code library failed to load.", Toast.LENGTH_SHORT).show();
//            status = false;
//        }
//
//        return status;
//    }

//    private native void initFilePath(String external_dir);
}