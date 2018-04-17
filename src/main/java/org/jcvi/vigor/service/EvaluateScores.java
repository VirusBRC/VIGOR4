package org.jcvi.vigor.service;

import java.util.Map;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.stereotype.Service;

@Service
public class EvaluateScores implements EvaluateModel {
	private float exonerateScoreFactor=1;
	private float startScoreFactor=1;
	private float splicingScoreFactor=1;
	private float stopScoreFactor=1;
	private float leakyStopScoreFactor=1;

	@Override
	public Model evaluate(Model model,VigorForm form) {
		
		Map<String,Double> scores = model.getScores();
		VigorConfiguration configuration = form.getConfiguration();
		if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.ScoreFactorExonerate))) {
			exonerateScoreFactor = Integer.parseInt(configuration.get(ConfigurationParameters.ScoreFactorExonerate));
		}
		if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.ScoreFactorStart))) {
			startScoreFactor = Integer.parseInt(configuration.get(ConfigurationParameters.ScoreFactorStart));
		}
		if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.ScoreFactorSplicing))) {
			splicingScoreFactor = Integer.parseInt(configuration.get(ConfigurationParameters.ScoreFactorSplicing));
		}
		if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.ScoreFactorStop))) {
			stopScoreFactor = Integer.parseInt(configuration.get(ConfigurationParameters.ScoreFactorStop));
		}
		if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.ScoreFactorLeakyStop))) {
			leakyStopScoreFactor = Integer.parseInt(configuration.get(ConfigurationParameters.ScoreFactorLeakyStop));
		}
		double exonerateScore=0;
		double startCodonScore=0;
		double splicingScore=0;
		double leakyStopScore=0;
		double stopScore=0;
		double totalScore=0;
		if(scores.get("exonerateScore")!=null){
		exonerateScore = scores.get("exonerateScore")*exonerateScoreFactor;
		}
		if(scores.get("startCodonScore")!=null){
	    startCodonScore = scores.get("startCodonScore")*startScoreFactor;
		}
		if(scores.get("leakyStopScore") !=null){
	    leakyStopScore = scores.get("leakyStopScore")*leakyStopScoreFactor;
		}
		if(scores.get("spliceScore")!=null){
	    splicingScore = scores.get("spliceScore")*splicingScoreFactor;
		}
		if(scores.get("stopCodonScore")!=null){
	    stopScore = scores.get("stopCodonScore")*stopScoreFactor;
		}		
		scores.put("exonerateScore",exonerateScore);
		scores.put("startCodonScore",startCodonScore);
		scores.put("leakyStopScore",leakyStopScore);
		scores.put("spliceScore",splicingScore);
		scores.put("stopCodonScore",stopScore);
		totalScore = exonerateScore + startCodonScore + leakyStopScore + splicingScore + stopScore;
		scores.put("totalScore", totalScore);
		return model;
	}
	
}
