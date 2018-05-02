package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.Appender;
import org.apache.logging.log4j.core.Layout;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.config.AppenderRef;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.config.LoggerConfig;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.jcvi.vigor.exception.VigorException;
import org.springframework.core.io.ClassPathResource;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.regex.Pattern;


public class VigorUtils {

	private static final Logger LOGGER = LogManager.getLogger ( VigorUtils.class );
	private static String BLAST_CLASSPATH = "vigorResources" + File.separator + "blast" + File.separator;
	private static String WINDOWS = "windows";
	private static String LINUX = "linux";

	private static Pattern hypPattern = Pattern.compile("^HYP\\b");

	public static boolean isNullOrEmpty(Object value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase(""))
			return false;
		return true;
	}

	public static String removeSpaces(String value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase("")) {
			return value.replaceAll("(?m)^[ \t]*\r?\n", "").trim();
		}
		return "";
	}

	/*public static Boolean isDefLine(String value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase("")) {
			return Pattern.compile("^>.*").matcher(value.trim()).matches();
		}
		return false;
	}*/

	/*public static Boolean isNucleotide(String value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase("")) {
			return Pattern.compile("^[a-zA-Z]*$").matcher(value.trim()).matches();
		}
		return false;
	}

	*/
	public static String getVigorWorkSpace() throws VigorException {
		String workingDirectory = null;
		for (String dirName: new String[] {
				System.getenv("VIGOR_TEMPDIR"),
				System.getProperty("vigor.tempdir"),
				"VigorWorkSpace"}) {

			if (dirName == null) {
				continue;
			}
			File theDir = new File(dirName);
			try {
				if (theDir.exists()) {
					if ( theDir.isDirectory() && theDir.canWrite() && theDir.canRead()) {
						workingDirectory = theDir.getPath();
					} else {
						LOGGER.error("{} is not a directory or is not readable and writable", theDir);
					}
				} else {
					if (theDir.mkdirs() ) {
						workingDirectory = theDir.getPath();
					} else {
						LOGGER.error("unable to create directory {}", theDir);
					}
				}
			} catch (SecurityException e) {
				LOGGER.error(e.getMessage(),e);
				throw new VigorException(String.format("unable to create working directory %s", theDir),e);
			}
			break;
		}
		if (workingDirectory == null) {
			throw new VigorException("unable to determine/create working directory");
		}
		return workingDirectory;
	}

	public static String getBlastCommand(String blastFilePath, String inputFilePath, String reference_db,
			String outputFilePath, String logFilePath) {
		return blastFilePath + " -p blastx -i " + inputFilePath + " -d " + reference_db
				+ " -e 1e-5 -M BLOSUM45 -g F -F \"\" -z 3000000 -v 0 -b 100 -m 7";
	}

	/*public static String getVirusDatabasePath() {
		String filepath="";
		try {
			File folder = new ClassPathResource ( "vigorResources/data3" ).getFile ();
			filepath=folder.getAbsolutePath ();

		}
		catch(IOException e){
			LOGGER.error ( e.getMessage (),e );
		}
		return filepath;
	}*/

	public static String getVigorParametersPath() {
		return Paths.get("vigorResources", "config", "vigor.ini").toString();
	}
	/*public static String getVirusSpecificParametersPath(){
        String virusParamsPath = "vigorResources" + File.separator + "config" + File.separator + "virusSpecificParams";
        return virusParamsPath;
    }*/

	public static String getBlastFilePath() throws IllegalArgumentException {
		if (OSValidator.isWindows()) {
			if (OSValidator.is64Bit()) {
				return BLAST_CLASSPATH + WINDOWS + File.separator + "64" + File.separator + "blastall";
			} else {
				return BLAST_CLASSPATH + WINDOWS + File.separator + "32" + File.separator + "blastall";
			}
		} else if (OSValidator.isUnix()) {
			if (OSValidator.is64Bit()) {
				return BLAST_CLASSPATH + LINUX + File.separator + "64" + File.separator + "blastall";
			} else {
				return BLAST_CLASSPATH + LINUX + File.separator + "32" + File.separator + "blastall";
			}
		} else if (OSValidator.isMac()) {
			/*
			 * if (OSValidator.is64Bit()) { return
			 * "blast"+File.pathSeparator+"64"+File.separator+"blastall"; } else
			 * { return
			 * "blast"+File.pathSeparator+"64"+File.separator+"blastall"; }
			 */
			throw new IllegalArgumentException("No Blast File Found for Mac OS");
		} else if (OSValidator.isSolaris()) {
			/*
			 * if (OSValidator.is64Bit()) { return
			 * "blast"+File.pathSeparator+"64"+File.pathSeparator+"blastall"; }
			 * else { return
			 * "blast"+File.pathSeparator+"64"+File.pathSeparator+"blastall"; }
			 */
			throw new IllegalArgumentException("No Blast File Found for Solaris OS");
		} else {
			System.out.println("Your OS is not support!!");
			throw new IllegalArgumentException("Your OS is not support!!");
		}
	}
	public static boolean is_Integer(String value){
		try {
			Integer.parseInt ( value );
			return true;
		}
		catch(NumberFormatException e)
		{
			return false;
		}
	}

	public static String removeQuotes(String in) {
		return in.replaceAll("^\"|\"$","");
	}

	public static String nameToLocus(String name, String prefix, boolean isPseudoGene) {
		if (prefix == null || prefix.isEmpty()) {
			return null;
		}
		String symbol = name;
		if (isPseudoGene || hypPattern.matcher(symbol).find()) {
			symbol = symbol.replace("-","p");
		} else {
			symbol = symbol.replace("-[0-9].*$","");
		}
		symbol = symbol.replace("[^a-zA-Z0-9-]","");
		return prefix + symbol;
	}

}
