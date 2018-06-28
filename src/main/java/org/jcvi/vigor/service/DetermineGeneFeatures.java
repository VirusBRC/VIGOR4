package org.jcvi.vigor.service;

import java.util.List;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;

public interface DetermineGeneFeatures {

    List<Model> determine ( Model model, VigorForm form ) throws ServiceException;
}
