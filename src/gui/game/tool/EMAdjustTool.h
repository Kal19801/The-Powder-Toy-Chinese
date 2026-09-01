#pragma once
#include "Tool.h"
#include "common/String.h"
#include <optional>

class GameModel;
class Simulation;

// Unified EM adjust tool: the six separate applet adjust modes (conductivity /
// permeability / current / dielectric / magnet direction / magnet strength)
// merged into one PROP-like drawing tool. A configuration window selects the
// property, the target value and one of three application modes:
//  - set:      the values under the brush become the target value
//  - add:      the target value is added, repeating every 0.2s while held
//  - subtract: the target value is subtracted, repeating every 0.2s while held
// The values are stored as EM cell overrides inside EMField, which only apply
// while the underlying particles still provide a matching material.
class EMAdjustTool: public Tool
{
public:
        struct Configuration
        {
                int property;    // one of EmAdjustProperty
                int applyMode;   // one of EmAdjustApply
                float value;     // target value; direction in degrees for MAG_DIR
                String valueStr;
        };

private:
        void Apply(Simulation *sim, ui::Point position);
        void ApplyBrush(Simulation *sim, Brush const &brush, ui::Point position);
        void SetConfiguration(std::optional<Configuration> newConfiguration);

        GameModel &gameModel;
        std::optional<Configuration> configuration;
        // frame counter for the 0.2s repeat rate of the add / subtract modes
        int repeatCounter = 0;
        // remembers whether the last add / subtract was applied, so that every
        // new stroke starts with an immediate application
        bool strokeActive = false;

        friend class EMAdjustWindow;

public:
        EMAdjustTool(GameModel &newGameModel):
                Tool(0, "EMADJ", ByteString("EM参数调整:以PROP笔的方式调整电磁场参数(电导/磁导/电流/介电/磁向/磁强)").FromUtf8(),
                        0x80BFFF_rgb, "DEFAULT_UI_EMADJUST", nullptr
                ), gameModel(newGameModel)
        {
                MenuSection = SC_EM;
        }

        void OpenWindow(Simulation *sim);
        void Click(Simulation * sim, Brush const &brush, ui::Point position) override { }
        void Draw(Simulation * sim, Brush const &brush, ui::Point position) override;
        void DrawLine(Simulation * sim, Brush const &brush, ui::Point position1, ui::Point position2, bool dragging) override;
        void DrawRect(Simulation * sim, Brush const &brush, ui::Point position1, ui::Point position2) override;
        void DrawFill(Simulation * sim, Brush const &brush, ui::Point position) override;

        std::optional<Configuration> GetConfiguration() const
        {
                return configuration;
        }

        void Select(int toolSelection) final override;
};
