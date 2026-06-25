package com.bitpop

import android.content.ContentValues
import android.net.Uri
import android.os.Build
import android.os.Bundle
import android.os.Environment
import android.provider.MediaStore
import android.provider.OpenableColumns
import android.util.Log
import android.view.View
import android.widget.TextView
import android.widget.Button
import android.widget.ProgressBar
import android.widget.Spinner
import android.widget.EditText
import android.widget.LinearLayout
import androidx.activity.result.contract.ActivityResultContracts
import androidx.appcompat.app.AppCompatActivity
import java.io.File

class MainActivity : AppCompatActivity() {

    external fun mapReads(indexPath: String, readsPath: String, outputPath: String): String
    external fun fastConMap(indexPaths: String, readsPath: String, outputPath: String, chunkPct: Int, chunkMin: Int, chunkMax: Int): String
    external fun buildIndex(fastaPath: String, outputPath: String, kmerSize: Int): String
    external fun getGenomeNames(indexPath: String): String

    companion object {
        init {
            System.loadLibrary("bit_pop")
        }
    }

    private var index1Path: String? = null
    private var index2Path: String? = null
    private var index3Path: String? = null
    private var readsPath: String? = null
    private var mode: String = "base"

    private val pickFile = registerForActivityResult(ActivityResultContracts.GetContent()) { uri: Uri? ->
        uri?.let {
            it.toFilePath()?.let { path: String ->
                currentPicker?.invoke(path)
            }
        }
    }

    private var currentPicker: ((String) -> Unit)? = null

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        val modeSpinner = findViewById<Spinner>(R.id.modeSpinner)
        val multiIndexLayout = findViewById<LinearLayout>(R.id.multiIndexLayout)
        val chunkingLayout = findViewById<LinearLayout>(R.id.chunkingLayout)
        val statusText = findViewById<TextView>(R.id.statusText)
        val progress = findViewById<ProgressBar>(R.id.progressBar)

        val modes = arrayOf("Base (single index)", "Fast-Con (multi-index)")
        val adapter = android.widget.ArrayAdapter(this, android.R.layout.simple_spinner_dropdown_item, modes)
        modeSpinner.adapter = adapter

        modeSpinner.onItemSelectedListener = object : android.widget.AdapterView.OnItemSelectedListener {
            override fun onItemSelected(p0: android.widget.AdapterView<*>?, p1: View?, position: Int, p3: Long) {
                mode = if (position == 0) "base" else "concon"
                multiIndexLayout.visibility = if (mode == "concon") LinearLayout.VISIBLE else LinearLayout.GONE
                chunkingLayout.visibility = if (mode == "concon") LinearLayout.VISIBLE else LinearLayout.GONE
            }
            override fun onNothingSelected(p0: android.widget.AdapterView<*>?) {}
        }

        findViewById<Button>(R.id.selectIndex1Button).setOnClickListener {
            currentPicker = { index1Path = it; findViewById<TextView>(R.id.index1PathText).text = it.substringAfterLast('/') }
            pickFile.launch("*/*")
        }

        findViewById<Button>(R.id.selectIndex2Button).setOnClickListener {
            currentPicker = { index2Path = it; findViewById<TextView>(R.id.index2PathText).text = it.substringAfterLast('/') }
            pickFile.launch("*/*")
        }

        findViewById<Button>(R.id.selectIndex3Button).setOnClickListener {
            currentPicker = { index3Path = it; findViewById<TextView>(R.id.index3PathText).text = it.substringAfterLast('/') }
            pickFile.launch("*/*")
        }

        findViewById<Button>(R.id.selectReadsButton).setOnClickListener {
            currentPicker = { readsPath = it; findViewById<TextView>(R.id.readsPathText).text = it.substringAfterLast('/') }
            pickFile.launch("*/*")
        }

        findViewById<Button>(R.id.runButton).setOnClickListener {
            val indexPath1 = index1Path
            val reads = readsPath

            if (indexPath1 == null) {
                statusText.text = "Select index 1 and reads"
                return@setOnClickListener
            }
            if (reads == null) {
                statusText.text = "Select reads file"
                return@setOnClickListener
            }

            progress.visibility = ProgressBar.VISIBLE
            statusText.text = "Running..."

            Thread {
                try {
                    val outputPath = File(cacheDir, "results.tsv").absolutePath

                    val result = if (mode == "base") {
                        mapReads(indexPath1!!, reads!!, outputPath)
                    } else {
                        val paths = mutableListOf(indexPath1!!)
                        index2Path?.let { paths.add(it) }
                        index3Path?.let { paths.add(it) }

                        val chunkPctEdit = findViewById<EditText>(R.id.chunkPctEdit)
                        val chunkPct = try { chunkPctEdit.text.toString().toInt() } catch (e: Exception) { 3 }

                        fastConMap(paths.joinToString(","), reads!!, outputPath, chunkPct, 20, 500)
                    }

                    val resultsFile = File(outputPath)
                    val preview = if (resultsFile.exists()) {
                        val lines = resultsFile.readLines().take(10)
                        lines.joinToString("\n")
                    } else {
                        "(empty)"
                    }

                    val downloadPath = saveToDownload(resultsFile, "bitpop_results.tsv")

                    runOnUiThread {
                        statusText.text = "$result\n\nSaved: $downloadPath\n\n--- Results (first 10) ---\n$preview"
                        progress.visibility = ProgressBar.GONE
                    }
                } catch (e: Exception) {
                    Log.e("BitPop", "Error", e)
                    runOnUiThread {
                        statusText.text = "Error: " + e.message
                        progress.visibility = ProgressBar.GONE
                    }
                }
            }.start()
        }
    }

    fun saveToDownload(source: File, filename: String): String {
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.Q) {
            val contentValues = ContentValues().apply {
                put(MediaStore.MediaColumns.DISPLAY_NAME, filename)
                put(MediaStore.MediaColumns.MIME_TYPE, "text/tab-separated-values")
                put(MediaStore.MediaColumns.RELATIVE_PATH, Environment.DIRECTORY_DOWNLOADS)
            }
            val uri = contentResolver.insert(MediaStore.Downloads.EXTERNAL_CONTENT_URI, contentValues)
            uri?.let {
                contentResolver.openOutputStream(it)?.use { output ->
                    source.inputStream().use { it -> it.copyTo(output) }
                }
                return "Download/$filename"
            }
        } else {
            val dest = File(Environment.getExternalStorageDirectory(), filename)
            source.copyTo(dest, overwrite = true)
            return dest.absolutePath
        }
        return "(save failed)"
    }

    fun Uri.toFilePath(): String? {
        return when (this.scheme) {
            "file" -> this.path
            "content" -> {
                val cursor = contentResolver.query(this, null, null, null, null)
                cursor?.use {
                    val nameIndex = it.getColumnIndex(OpenableColumns.DISPLAY_NAME)
                    if (it.moveToFirst() && nameIndex != -1) {
                        val fileName = it.getString(nameIndex)
                        val destFile = File(cacheDir, fileName)
                        contentResolver.openInputStream(this)?.use { input ->
                            destFile.outputStream().use { output ->
                                input.copyTo(output)
                            }
                        }
                        destFile.absolutePath
                    } else null
                }
            }
            else -> null
        }
    }
}