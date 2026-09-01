#include "EMAdjustTool.h"
#include "prefs/GlobalPrefs.h"
#include "gui/Style.h"
#include "gui/game/Brush.h"
#include "gui/game/GameModel.h"
#include "gui/interface/Window.h"
#include "gui/interface/Button.h"
#include "gui/interface/DropDown.h"
#include "gui/interface/Textbox.h"
#include "simulation/Simulation.h"
#include "simulation/EMField.h"
#include "simulation/CoordStack.h"
#include "graphics/Graphics.h"
#include "Format.h"
#include <SDL.h>
#include <cmath>

// Unified EM adjust tool window, laid out like the PROP tool window: a property
// dropdown, a value textbox, plus a mode dropdown (set / add / subtract).
class EMAdjustWindow: public ui::Window
{
        void Update();

public:
        ui::DropDown * property;
        ui::DropDown * applyMode;
        ui::Textbox * textField;
        ui::Button * okayButton;
        EMAdjustTool * tool;
        std::optional<EMAdjustTool::Configuration> configuration;

        EMAdjustWindow(EMAdjustTool *tool_);
        void OnTryExit(ExitMethod method) override;
        void OnDraw() override;
        void OnKeyPress(int key, int scan, bool repeat, bool shift, bool ctrl, bool alt) override;
        virtual ~EMAdjustWindow() {}
};

EMAdjustWindow::EMAdjustWindow(EMAdjustTool *tool_):
        ui::Window(ui::Point(-1, -1), ui::Point(200, 107)),
        tool(tool_)
{
        ui::Label * messageLabel = new ui::Label(ui::Point(4, 5), ui::Point(Size.X-8, 14), ByteString("EM参数调整工具").FromUtf8());
        messageLabel->SetTextColour(style::Colour::InformationTitle);
        messageLabel->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
        messageLabel->Appearance.VerticalAlign = ui::Appearance::AlignTop;
        AddComponent(messageLabel);

        ui::Button * okayButton = new ui::Button(ui::Point(0, Size.Y-17), ui::Point(Size.X, 17), "OK");
        okayButton->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
        okayButton->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
        okayButton->Appearance.BorderInactive = ui::Colour(200, 200, 200);
        okayButton->SetActionCallback({ [this] {
                CloseActiveWindow();
                if (configuration)
                {
                        tool->SetConfiguration(configuration);
                }
                SelfDestruct();
        } });
        AddComponent(okayButton);
        SetOkayButton(okayButton);
        this->okayButton = okayButton;

        property = new ui::DropDown(ui::Point(8, 25), ui::Point(Size.X-16, 16));
        property->SetActionCallback({ [this] {
                Update();
        } });
        AddComponent(property);
        property->AddOption({ ByteString("电导率 Conductivity (0~1)").FromUtf8(), EMADJP_CONDUCT });
        property->AddOption({ ByteString("磁导率 Permeability (0.05~32)").FromUtf8(), EMADJP_PERM });
        property->AddOption({ ByteString("电流 Current (-2~2)").FromUtf8(), EMADJP_J });
        property->AddOption({ ByteString("介电常数 Dielectric (1~191)").FromUtf8(), EMADJP_MEDIUM });
        property->AddOption({ ByteString("磁方向 Mag Dir (0~360度)").FromUtf8(), EMADJP_MAG_DIR });
        property->AddOption({ ByteString("磁强度 Mag Strength (0~2)").FromUtf8(), EMADJP_MAG_STR });

        applyMode = new ui::DropDown(ui::Point(8, 45), ui::Point(Size.X-16, 16));
        applyMode->SetActionCallback({ [this] {
                Update();
        } });
        AddComponent(applyMode);
        applyMode->AddOption({ ByteString("设置模式:设为目标值").FromUtf8(), EMADJA_SET });
        applyMode->AddOption({ ByteString("加法模式:加上目标值(0.2s/次)").FromUtf8(), EMADJA_ADD });
        applyMode->AddOption({ ByteString("减法模式:减去目标值(0.2s/次)").FromUtf8(), EMADJA_SUB });

        textField = new ui::Textbox(ui::Point(8, 65), ui::Point(Size.X-16, 16), "", ByteString("[值]").FromUtf8());
        textField->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
        textField->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
        textField->SetActionCallback({ [this]() {
                Update();
        } });
        AddComponent(textField);

        auto &prefs = GlobalPrefs::Ref();
        int propertyIndex = prefs.Get("EMAdj.Property", int(EMADJP_CONDUCT));
        int modeIndex = prefs.Get("EMAdj.Mode", int(EMADJA_SET));
        auto valueString = prefs.Get("EMAdj.Value", String(""));
        if (propertyIndex < 0 || propertyIndex >= EMADJP_COUNT)
        {
                propertyIndex = EMADJP_CONDUCT;
        }
        if (modeIndex < 0 || modeIndex >= EMADJA_SUB)
        {
                modeIndex = EMADJA_SET;
        }
        property->SetOption(propertyIndex);
        applyMode->SetOption(modeIndex);
        textField->SetText(valueString);

        FocusComponent(textField);
        Update();

        MakeActiveWindow();
}

void EMAdjustWindow::Update()
{
        configuration.reset();
        auto valueStr = textField->GetText();
        float value = 0;
        bool ok = true;
        try
        {
                value = valueStr.ToNumber<float>();
        }
        catch (const std::exception &)
        {
                ok = false;
        }
        // magnet direction is given in degrees and may be any number
        if (ok || property->GetOption().second == EMADJP_MAG_DIR)
        {
                configuration = EMAdjustTool::Configuration{
                        property->GetOption().second,
                        applyMode->GetOption().second,
                        ok ? value : 0.0f,
                        valueStr,
                };
        }
        auto haveConfiguration = bool(configuration);
        okayButton->Enabled = haveConfiguration;
        textField->SetTextColour(haveConfiguration ? ui::Colour(255, 255, 255) : style::Colour::ErrorTitle);
}

void EMAdjustWindow::OnTryExit(ExitMethod method)
{
        CloseActiveWindow();
        SelfDestruct();
}

void EMAdjustWindow::OnDraw()
{
        Graphics * g = GetGraphics();

        g->DrawFilledRect(RectSized(Position - Vec2{ 1, 1 }, Size + Vec2{ 2, 2 }), 0x000000_rgb);
        g->DrawRect(RectSized(Position, Size), 0xC8C8C8_rgb);
}

void EMAdjustWindow::OnKeyPress(int key, int scan, bool repeat, bool shift, bool ctrl, bool alt)
{
        if (key == SDLK_UP)
        {
                property->SetOption(property->GetOption().second-1);
                Update();
        }
        else if (key == SDLK_DOWN)
        {
                property->SetOption(property->GetOption().second+1);
                Update();
        }
}

void EMAdjustTool::SetConfiguration(std::optional<Configuration> newConfiguration)
{
        configuration = newConfiguration;
        if (configuration)
        {
                auto &prefs = GlobalPrefs::Ref();
                Prefs::DeferWrite dw(prefs);
                prefs.Set("EMAdj.Property", configuration->property);
                prefs.Set("EMAdj.Mode", configuration->applyMode);
                prefs.Set("EMAdj.Value", configuration->valueStr);
        }
}

void EMAdjustTool::OpenWindow(Simulation *sim)
{
        new EMAdjustWindow(this);
}

// apply the configured property change to one EM cell
void EMAdjustTool::Apply(Simulation *sim, ui::Point position)
{
        if (position.X < 0 || position.X >= XRES || position.Y < 0 || position.Y >= YRES || !configuration)
        {
                return;
        }
        auto *emf = sim->GetEMField();
        if (!emf || !emf->enabled)
        {
                return;
        }
        int gi = emf->CellIndex(position.X, position.Y);
        if (emf->ApplyEMProperty(configuration->property, configuration->applyMode, gi, configuration->value))
        {
                // the wave equation needs to know about the material change right away
                emf->CalcBoundaries();
        }
}

// apply the configured property change to every cell under the brush, sharing
// the 0.2s repeat throttle of the add / subtract modes between all shapes
void EMAdjustTool::ApplyBrush(Simulation *sim, Brush const &cBrush, ui::Point position)
{
        if (!configuration)
        {
                return;
        }
        if (configuration->applyMode == EMADJA_SET)
        {
                // set mode applies continuously, exactly like the PROP tool
                for (ui::Point off : cBrush)
                {
                        ui::Point coords = position + off;
                        if (coords.X >= 0 && coords.Y >= 0 && coords.X < XRES && coords.Y < YRES)
                        {
                                Apply(sim, coords);
                        }
                }
        }
        else
        {
                // add / subtract repeat at 0.2s intervals (12 frames at 60 fps);
                // the first application of a stroke is immediate
                if (!strokeActive || repeatCounter <= 0)
                {
                        for (ui::Point off : cBrush)
                        {
                                ui::Point coords = position + off;
                                if (coords.X >= 0 && coords.Y >= 0 && coords.X < XRES && coords.Y < YRES)
                                {
                                        Apply(sim, coords);
                                }
                        }
                        strokeActive = true;
                        repeatCounter = 12;
                }
                else
                {
                        repeatCounter--;
                }
        }
}

void EMAdjustTool::Draw(Simulation *sim, Brush const &cBrush, ui::Point position)
{
        ApplyBrush(sim, cBrush, position);
}

void EMAdjustTool::DrawLine(Simulation *sim, Brush const &cBrush, ui::Point position, ui::Point position2, bool dragging)
{
        int x1 = position.X, y1 = position.Y, x2 = position2.X, y2 = position2.Y;
        bool reverseXY = abs(y2-y1) > abs(x2-x1);
        int x, y, dx, dy, sy, rx = cBrush.GetRadius().X, ry = cBrush.GetRadius().Y;
        float e = 0.0f, de;
        if (reverseXY)
        {
                y = x1;
                x1 = y1;
                y1 = y;
                y = x2;
                x2 = y2;
                y2 = y;
        }
        if (x1 > x2)
        {
                y = x1;
                x1 = x2;
                x2 = y;
                y = y1;
                y1 = y2;
                y2 = y;
        }
        dx = x2 - x1;
        dy = abs(y2 - y1);
        if (dx)
        {
                de = dy/(float)dx;
        }
        else
        {
                de = 0.0f;
        }
        y = y1;
        sy = (y1<y2) ? 1 : -1;
        for (x=x1; x<=x2; x++)
        {
                if (reverseXY)
                {
                        Draw(sim, cBrush, ui::Point(y, x));
                }
                else
                {
                        Draw(sim, cBrush, ui::Point(x, y));
                }
                e += de;
                if (e >= 0.5f)
                {
                        y += sy;
                        if (!(rx+ry) && ((y1<y2) ? (y<=y2) : (y>=y2)))
                        {
                                if (reverseXY)
                                {
                                        Draw(sim, cBrush, ui::Point(y, x));
                                }
                                else
                                {
                                        Draw(sim, cBrush, ui::Point(x, y));
                                }
                        }
                        e -= 1.0f;
                }
        }
}

void EMAdjustTool::DrawRect(Simulation *sim, Brush const &cBrush, ui::Point position, ui::Point position2)
{
        int x1 = position.X, y1 = position.Y, x2 = position2.X, y2 = position2.Y;
        int i, j;
        if (x1>x2)
        {
                i = x2;
                x2 = x1;
                x1 = i;
        }
        if (y1>y2)
        {
                j = y2;
                y2 = y1;
                y1 = j;
        }
        for (j=y1; j<=y2; j++)
        {
                for (i=x1; i<=x2; i++)
                {
                        ApplyBrush(sim, cBrush, ui::Point(i, j));
                }
        }
}

// flood fill over EM cells: spreads the property change through cells of the
// same material type as the starting cell (the EM cell pendant of flood_prop)
void EMAdjustTool::DrawFill(Simulation *sim, Brush const &cBrush, ui::Point position)
{
        if (!configuration || position.X < 0 || position.X >= XRES || position.Y < 0 || position.Y >= YRES)
        {
                return;
        }
        auto *emf = sim->GetEMField();
        if (!emf || !emf->enabled)
        {
                return;
        }
        int start = emf->CellIndex(position.X, position.Y);
        int startType = EMField::CellTypeOf(emf->cells[start]);
        if (startType == EMCT_NONE)
        {
                return;
        }
        CoordStack stack;
        stack.push(position.X, position.Y);
        int cellsChanged = 0;
        while (stack.getSize() && cellsChanged < emf->gw * emf->gh)
        {
                int x, y;
                stack.pop(x, y);
                if (x < 0 || y < 0 || x >= XRES || y >= YRES)
                {
                        continue;
                }
                int gi = emf->CellIndex(x, y);
                if (EMField::CellTypeOf(emf->cells[gi]) != startType)
                {
                        continue;
                }
                Apply(sim, ui::Point(x, y));
                cellsChanged++;
                stack.push(x+1, y);
                stack.push(x-1, y);
                stack.push(x, y+1);
                stack.push(x, y-1);
        }
}

void EMAdjustTool::Select(int toolSelection)
{
        OpenWindow(gameModel.GetSimulation());
}
