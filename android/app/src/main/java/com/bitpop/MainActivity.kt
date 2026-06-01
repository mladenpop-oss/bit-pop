package com.bitpop

import android.os.Bundle
import android.widget.TextView
import android.widget.Button
import android.widget.ProgressBar
import androidx.appcompat.app.AppCompatActivity
import java.io.File

class MainActivity : AppCompatActivity() {

    external fun mapReads(indexPath: String, readsPath: String, outputPath: String): String
    external fun buildIndex(fastaPath: String, outputPath: String, kmerSize: Int): String
    external fun getGenomeNames(indexPath: String): String

    companion object {
        init {
            System.loadLibrary("bitpop")
        }
    }

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        val statusText = findViewById<TextView>(R.id.statusText)
        val mapButton = findViewById<Button>(R.id.mapButton)
        val progress = findViewById<ProgressBar>(R.id.progressBar)

        val docsDir = filesDir.absolutePath
        statusText.text = "Bit-Pop Android\nData dir: $docsDir\n\nLoad index and reads files to begin."

        mapButton.setOnClickListener {
            progress.visibility = ProgressBar.VISIBLE
            statusText.text = "Mapping reads..."

            Thread {
                try {
                    val indexPath = "$docsDir/ebola_multi.bitpop"
                    val readsPath = "$docsDir/ebola_reads.fastq"
                    val outputPath = "$docsDir/results.tsv"

                    if (!File(indexPath).exists() || !File(readsPath).exists()) {
                        runOnUiThread {
                            statusText.text = "Error: Index or reads file not found.\nCopy files to: $docsDir"
                            progress.visibility = ProgressBar.GONE
                        }
                        return@Thread
                    }

                    val result = mapReads(indexPath, readsPath, outputPath)
                    val names = getGenomeNames(indexPath)

                    runOnUiThread {
                        statusText.text = "Result: $result\n\nGenomes: $names"
                        progress.visibility = ProgressBar.GONE
                    }
                } catch (e: Exception) {
                    runOnUiThread {
                        statusText.text = "Error: ${e.message}"
                        progress.visibility = ProgressBar.GONE
                    }
                }
            }.start()
        }
    }
}
