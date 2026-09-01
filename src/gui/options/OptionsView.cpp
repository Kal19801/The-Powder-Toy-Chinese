#include "OptionsView.h"
#include "Format.h"
#include "OptionsController.h"
#include "OptionsModel.h"
#include "common/clipboard/Clipboard.h"
#include "common/platform/Platform.h"
#include "graphics/Graphics.h"
#include "graphics/Renderer.h"
#include "gui/Style.h"
#include "simulation/ElementDefs.h"
#include "simulation/SimulationSettings.h"
#include "simulation/Air.h"
#include "client/Client.h"
#include "gui/credits/Credits.h"
#include "gui/dialogues/ConfirmPrompt.h"
#include "gui/dialogues/InformationMessage.h"
#include "gui/interface/Button.h"
#include "gui/interface/Checkbox.h"
#include "gui/interface/DropDown.h"
#include "gui/interface/Engine.h"
#include "gui/interface/Label.h"
#include "gui/interface/Separator.h"
#include "gui/interface/Textbox.h"
#include "gui/interface/DirectionSelector.h"
#include "PowderToySDL.h"
#include "Config.h"
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <SDL.h>

class DirectionSelector : public ui::Window
{
        std::function<void (float, float)> done;

        void OnTryExit(ExitMethod method) override
        {
                CloseActiveWindow();
                SelfDestruct();
        }

        void OnDraw() override
        {
                Graphics * g = GetGraphics();

                g->DrawFilledRect(RectSized(Position - Vec2{ 1, 1 }, Size + Vec2{ 2, 2 }), 0x000000_rgb);
                g->DrawRect(RectSized(Position, Size), 0xC8C8C8_rgb);
        }

        ui::DirectionSelector * direction;
        ui::Label * labelValues;

public:
        DirectionSelector(ui::Point position, float scale, int radius, float x, float y, String label, std::function<void (float, float)> newDone):
                ui::Window(position, ui::Point((radius * 5 / 2) + 20, (radius * 5 / 2) + 75)),
                done(newDone),
                direction(new ui::DirectionSelector(ui::Point(10, 32), scale, radius, radius / 4, 2, 5))
        {
                ui::Label * tempLabel = new ui::Label(ui::Point(4, 1), ui::Point(Size.X - 8, 22), label);
                tempLabel->SetTextColour(style::Colour::InformationTitle);
                tempLabel->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
                tempLabel->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                AddComponent(tempLabel);

                auto * tempSeparator = new ui::Separator(ui::Point(0, 22), ui::Point(Size.X, 1));
                AddComponent(tempSeparator);

                labelValues = new ui::Label(ui::Point(0, (radius * 5 / 2) + 37), ui::Point(Size.X, 16), String::Build(Format::Precision(1), "X:", x, " Y:", y, " Total:", std::hypot(x, y)));
                labelValues->Appearance.HorizontalAlign = ui::Appearance::AlignCentre;
                labelValues->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                AddComponent(labelValues);

                direction->SetValues(x, y);
                direction->SetUpdateCallback([this](float x, float y) {
                        labelValues->SetText(String::Build(Format::Precision(1), "X:", x, " Y:", y, " Total:", std::hypot(x, y)));
                });
                direction->SetSnapPoints(5, 5, 2);
                AddComponent(direction);

                ui::Button * okayButton = new ui::Button(ui::Point(0, Size.Y - 17), ui::Point(Size.X, 17), "OK");
                okayButton->Appearance.HorizontalAlign = ui::Appearance::AlignCentre;
                okayButton->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                okayButton->Appearance.BorderInactive = ui::Colour(200, 200, 200);
                okayButton->SetActionCallback({ [this] {
                        done(direction->GetXValue(), direction->GetYValue());
                        CloseActiveWindow();
                        SelfDestruct();
                } });
                AddComponent(okayButton);
                SetOkayButton(okayButton);

                MakeActiveWindow();
        }
};

OptionsView::OptionsView() : ui::Window(ui::Point(-1, -1), ui::Point(320, 340))
{
        auto autoWidth = [this](ui::Component *c, int extra) {
                c->Size.X = Size.X - c->Position.X - 12 - extra;
        };
        
        {
                auto *label = new ui::Label(ui::Point(4, 1), ui::Point(Size.X-8, 22), ByteString("设置").FromUtf8());
                label->SetTextColour(style::Colour::InformationTitle);
                label->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
                label->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                autoWidth(label, 0);
                AddComponent(label);
        }

        auto *tmpSeparator = new ui::Separator(ui::Point(0, 22), ui::Point(Size.X, 1));
        AddComponent(tmpSeparator);

        scrollPanel = new ui::ScrollPanel(ui::Point(1, 23), ui::Point(Size.X-2, Size.Y-39));
        
        AddComponent(scrollPanel);

        int currentY = 8;
        auto addLabel = [this, &currentY, &autoWidth](int indent, String text) {
                auto *label = new ui::Label(ui::Point(22 + indent * 15, currentY), ui::Point(1, 16), "");
                autoWidth(label, 0);
                label->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
                label->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                label->SetMultiline(true);
                label->SetText("\bg" + text); // stupid hack because autoWidth just changes Size.X and that doesn't update the text wrapper
                label->AutoHeight();
                scrollPanel->AddChild(label);
                currentY += label->Size.Y - 1;
                return label;
        };
        auto addCheckbox = [this, &currentY, &autoWidth, &addLabel](int indent, String text, String info, std::function<void ()> action) {
                auto *checkbox = new ui::Checkbox(ui::Point(8 + indent * 15, currentY), ui::Point(1, 16), text, "");
                autoWidth(checkbox, 0);
                checkbox->SetActionCallback({ action });
                currentY += 14;
                if (info.size())
                {
                        addLabel(indent, info);
                }
                currentY += 4;
                scrollPanel->AddChild(checkbox);
                return checkbox;
        };
        auto addDropDown = [this, &currentY, &autoWidth](String info, std::vector<std::pair<String, int>> options, std::function<void ()> action) {
                auto *dropDown = new ui::DropDown(ui::Point(Size.X - 95, currentY), ui::Point(80, 16));
                scrollPanel->AddChild(dropDown);
                for (auto &option : options)
                {
                        dropDown->AddOption(option);
                }
                dropDown->SetActionCallback({ action });
                auto *label = new ui::Label(ui::Point(8, currentY), ui::Point(Size.X - 96, 16), info);
                label->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
                label->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                scrollPanel->AddChild(label);
                autoWidth(label, 85);
                currentY += 20;
                return dropDown;
        };
        auto addSeparator = [this, &currentY]() {
                currentY += 6;
                auto *separator = new ui::Separator(ui::Point(0, currentY), ui::Point(Size.X, 1));
                scrollPanel->AddChild(separator);
                currentY += 11;
        };
        auto addTextboxWithPreview = [this, &currentY](String info, bool addPreview, std::function<void (String value, bool defocus)> updateFunc) {
                auto *textbox = new ui::Textbox(ui::Point(Size.X-95, currentY), ui::Point(80, 16));
                textbox->SetActionCallback({ [textbox, updateFunc] {
                        updateFunc(textbox->GetText(), false);
                } });
                textbox->SetDefocusCallback({ [textbox, updateFunc] {
                        updateFunc(textbox->GetText(), true);
                } });
                textbox->SetLimit(9);
                scrollPanel->AddChild(textbox);
                ui::Button *preview{};
                if (addPreview)
                {
                        textbox->Size.X -= 20;
                        preview = new ui::Button(ui::Point(Size.X-31, currentY), ui::Point(16, 16), "", "Preview");
                        scrollPanel->AddChild(preview);
                }
                auto *label = new ui::Label(ui::Point(8, currentY), ui::Point(Size.X-105, 16), info);
                label->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
                label->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                scrollPanel->AddChild(label);
                currentY += 20;
                return std::make_pair(textbox, preview);
        };

        heatSimulation = addCheckbox(0, ByteString("热模拟 \bg34.0版本后加入").FromUtf8(), ByteString("\bg 关闭此选项可能导致一些奇怪的问题").FromUtf8(), [this] {
                c->SetHeatSimulation(heatSimulation->GetChecked());
        });
        newtonianGravity = addCheckbox(0, ByteString("牛顿引力模拟 \bg48.0版本后加入").FromUtf8(), ByteString("\bg 可能会降低游戏运行的效能").FromUtf8(), [this] {
                c->SetNewtonianGravity(newtonianGravity->GetChecked());
        });
        ambientHeatSimulation = addCheckbox(0, ByteString("环境热模拟 \bg50.0版本后加入").FromUtf8(), ByteString("\bg 关闭此项可能导致一些沙盘不能正常运行").FromUtf8(), [this] {
                c->SetAmbientHeatSimulation(ambientHeatSimulation->GetChecked());
        });
        waterEqualisation = addCheckbox(0, ByteString("连通器模拟 \bg61.0版本后加入 ").FromUtf8(), ByteString("\bg 有大量液体存在时会降低游戏运行的效能").FromUtf8(), [this] {
                c->SetWaterEqualisation(waterEqualisation->GetChecked());
        });
        airMode = addDropDown(ByteString("空气模拟模式").FromUtf8(), {
                { ByteString("开启").FromUtf8(), AIR_ON },
                { ByteString("关闭压力").FromUtf8(), AIR_PRESSUREOFF },
                { ByteString("关闭速度").FromUtf8(), AIR_VELOCITYOFF },
                { ByteString("关闭").FromUtf8(), AIR_OFF },
                { ByteString("更新停止").FromUtf8(), AIR_NOUPDATE },
        }, [this] {
                c->SetAirMode(airMode->GetOption().second);
        });
        std::tie(ambientAirTemp, ambientAirTempPreview) = addTextboxWithPreview("Ambient air temperature", true, [this](String value, bool defocus) {
                UpdateAirTemp(value, defocus);
        });
        std::tie(edgePressure, edgePressurePreview) = addTextboxWithPreview("Ambient air pressure", true, [this](String value, bool defocus) {
                UpdateEdgePressure(value, defocus);
        });
        {
                edgeVelocityChange = new ui::Button(ui::Point(Size.X-95, currentY), ui::Point(80, 16), "Change");
                scrollPanel->AddChild(edgeVelocityChange);
                edgeVelocityChange->SetActionCallback({ [this] {
                        new DirectionSelector(ui::Point(-1, -1), 0.05f, 40, edgeVelocityX, edgeVelocityY, "Ambient air velocity", [this](float x, float y) {
                                c->SetEdgeVelocityX(x);
                                c->SetEdgeVelocityY(y);
                        });
                } });
                auto *label = new ui::Label(ui::Point(8, currentY), ui::Point(Size.X-96, 16), "Ambient air velocity");
                label->Appearance.HorizontalAlign = ui::Appearance::AlignLeft;
                label->Appearance.VerticalAlign = ui::Appearance::AlignMiddle;
                scrollPanel->AddChild(label);
                currentY+=20;
        }
        vorticityCoeff = addTextboxWithPreview("Vorticity confinement", false, [this](String value, bool defocus) {
                UpdateVorticityCoeff(value, defocus);
        }).first;
        convectionMode = addDropDown("Air heat convection mode", {
                { "None", AIRC_NONE },
                { "Legacy", AIRC_LEGACY },
                { "Boussinesq", AIRC_BOUSSINESQ },
        }, [this] {
                c->SetConvectionMode(convectionMode->GetOption().second);
        });
        gravityMode = addDropDown("Gravity simulation mode", {
                { "Vertical", GRAV_VERTICAL },
                { "Off", GRAV_OFF },
                { "Radial", GRAV_RADIAL },
                { "Custom", GRAV_CUSTOM },
        }, [this] {
                c->SetGravityMode(gravityMode->GetOption().second);
                if (gravityMode->GetOption().second == 3)
                {
                        new DirectionSelector(ui::Point(-1, -1), 0.05f, 40, customGravityX, customGravityY, "Custom Gravity", [this](float x, float y) {
                                c->SetCustomGravityX(x);
                                c->SetCustomGravityY(y);
                        });
                }
        });
        edgeMode = addDropDown(ByteString("边界模式").FromUtf8(), {
                { ByteString("虚空").FromUtf8(), EDGE_VOID },
                { ByteString("固体").FromUtf8(), EDGE_SOLID },
                { ByteString("循环").FromUtf8(), EDGE_LOOP },
        }, [this] {
                c->SetEdgeMode(edgeMode->GetOption().second);
        });
        temperatureScale = addDropDown(ByteString("温度单位").FromUtf8(), {
                { ByteString("开尔文").FromUtf8(), TEMPSCALE_KELVIN },
                { ByteString("摄氏度").FromUtf8(), TEMPSCALE_CELSIUS },
                { ByteString("华氏度").FromUtf8(), TEMPSCALE_FAHRENHEIT },
        }, [this] {
                c->SetTemperatureScale(TempScale(temperatureScale->GetOption().second));
        });
        if (FORCE_WINDOW_FRAME_OPS != forceWindowFrameOpsHandheld)
        {
                addSeparator();
                language = addDropDown(ByteString("\bg语言").FromUtf8(), {  
                { ByteString("简体中文").FromUtf8(), 0 },  
                { ByteString("English").FromUtf8(), 1 },  
                }, [this] {  
                c->SetLanguage(language->GetOption().second);  
                });  
  
                std::vector<std::pair<String, int>> options;
                int currentScale = ui::Engine::Ref().GetScale();
                int scaleIndex = 1;
                bool currentScaleValid = false;
                do
                {
                        if (currentScale == scaleIndex)
                        {
                                currentScaleValid = true;
                        }
                        options.push_back({ String::Build(scaleIndex), scaleIndex });
                        scaleIndex += 1;
                }
                while (desktopWidth >= GetGraphics()->Size().X * scaleIndex && desktopHeight >= GetGraphics()->Size().Y * scaleIndex);
                if (!currentScaleValid)
                {
                        options.push_back({ ByteString("当前").FromUtf8(), currentScale });
                }
                scale = addDropDown(ByteString("\bg缩放屏幕的窗口比例因子").FromUtf8(), options, [this] {
                        c->SetScale(scale->GetOption().second);
                });
        }
        if (FORCE_WINDOW_FRAME_OPS == forceWindowFrameOpsNone)
        {
                resizable = addCheckbox(0, ByteString("可调整大小 \bg - 允许调整大小和最大化窗口").FromUtf8(), "", [this] {
                        c->SetResizable(resizable->GetChecked());
                });
                fullscreen = addCheckbox(0, ByteString("全屏 \bg - 进入全屏模式").FromUtf8(), "", [this] {
                        c->SetFullscreen(fullscreen->GetChecked());
                });
                changeResolution = addCheckbox(1, ByteString("设置最佳屏幕分辨率").FromUtf8(), "", [this] {
                        c->SetChangeResolution(changeResolution->GetChecked());
                });
                forceIntegerScaling = addCheckbox(1, ByteString("强制整数倍缩放 \bg - 不再那么模糊").FromUtf8(), "", [this] {
                        c->SetForceIntegerScaling(forceIntegerScaling->GetChecked());
                });
        }
        blurryScaling = addCheckbox(0, ByteString("模糊缩放 \bg - 启用此项加强在超大屏幕上效果").FromUtf8(), "", [this] {
                c->SetBlurryScaling(blurryScaling->GetChecked());
        });
        addSeparator();
        if (ALLOW_QUIT)
        {
                fastquit = addCheckbox(0, ByteString("快速退出").FromUtf8(), ByteString("点击关闭按钮时总是完全退出游戏").FromUtf8(), [this] {
                        c->SetFastQuit(fastquit->GetChecked());
                });
                globalQuit = addCheckbox(0, ByteString("全局退出快捷键").FromUtf8(), ByteString("使Ctrl+q 全局有效").FromUtf8(), [this] {
                        c->SetGlobalQuit(globalQuit->GetChecked());
                });
        }
        showAvatars = addCheckbox(0, ByteString("显示头像").FromUtf8(), ByteString("禁用此项可减少使用的网络带宽").FromUtf8(), [this] {
                c->SetShowAvatars(showAvatars->GetChecked());
        });
        momentumScroll = addCheckbox(0, ByteString("加速/旧版滚动").FromUtf8(), ByteString("启用此项时将步进滚动改为加速").FromUtf8(), [this] {
                c->SetMomentumScroll(momentumScroll->GetChecked());
        });
        mouseClickRequired = addCheckbox(0, ByteString("置顶类别").FromUtf8(), ByteString("启用此项时将滑动切换类别改为点击").FromUtf8(), [this] {
                c->SetMouseClickrequired(mouseClickRequired->GetChecked());
        });
        includePressure = addCheckbox(0, ByteString("压力数据").FromUtf8(), ByteString("沙盘,图章,复制时保存压力数据").FromUtf8(), [this] {
                c->SetIncludePressure(includePressure->GetChecked());
        });
        perfectCircle = addCheckbox(0, ByteString("完美的圆").FromUtf8(), ByteString("由Notch创造的最完美的圆").FromUtf8(), [this] {
                c->SetPerfectCircle(perfectCircle->GetChecked());
        });
        graveExitsConsole = addCheckbox(0, ByteString("Esc键下方的键用于退出控制台").FromUtf8(), ByteString("如果该键是0,请禁用此功能").FromUtf8(), [this] {
                c->SetGraveExitsConsole(graveExitsConsole->GetChecked());
        });
        if constexpr (PLATFORM_CLIPBOARD)
        {
                auto indent = 0;
                nativeClipoard = addCheckbox(indent, ByteString("使用全局剪贴板").FromUtf8(), ByteString("允许跨TPT实例进行复制和粘贴").FromUtf8(), [this] {
                        c->SetNativeClipoard(nativeClipoard->GetChecked());
                });
                currentY -= 4; // temporarily undo the currentY += 4 at the end of addCheckbox
                if (auto extra = Clipboard::Explanation())
                {
                        addLabel(indent, "\bg" + *extra);
                }
                currentY += 4; // and then undo the undo
        }
        threadedRendering = addCheckbox(0, ByteString("启用独立渲染线程").FromUtf8(), ByteString("使用特效时可能提升帧率").FromUtf8(), [this] {
                c->SetThreadedRendering(threadedRendering->GetChecked());
        });
        decoSpace = addDropDown(ByteString("\bg装饰工具使用的颜色空间").FromUtf8(), {
                { "sRGB", DECOSPACE_SRGB },
                { "Linear", DECOSPACE_LINEAR },
                { "Gamma 2.2", DECOSPACE_GAMMA22 },
                { "Gamma 1.8", DECOSPACE_GAMMA18 },
        }, [this] {
                c->SetDecoSpace(decoSpace->GetOption().second);
        });

        // TM-mode EM field section (EMWave2 port)
        addSeparator();
        emEnabled = addCheckbox(0, ByteString("电磁场模拟 \bg- TM模电磁波").FromUtf8(), ByteString("\bg 金属反射波, 玻璃折射, 火花与EMW辐射波, 波可加热/感应金属").FromUtf8(), [this] {
                c->SetEMEnabled(emEnabled->GetChecked());
        });
        emCellSize = addDropDown(ByteString("电磁场网格大小").FromUtf8(), {
                { ByteString("2 像素 (306x192)").FromUtf8(), 2 },
                { ByteString("4 像素 (153x96)").FromUtf8(), 4 },
                { ByteString("8 像素 (76x48)").FromUtf8(), 8 },
                { ByteString("16 像素 (38x24)").FromUtf8(), 16 },
        }, [this] {
                c->SetEMCellSize(emCellSize->GetOption().second);
        });
        emViewMode = addDropDown(ByteString("电磁场显示模式").FromUtf8(), {
                { ByteString("关闭").FromUtf8(), EMVIEW_OFF },
                { ByteString("电场 E").FromUtf8(), EMVIEW_E },
                { ByteString("磁场 B").FromUtf8(), EMVIEW_B },
                { ByteString("磁感线").FromUtf8(), EMVIEW_B_LINES },
                { ByteString("磁场强度 |B|").FromUtf8(), EMVIEW_B_STRENGTH },
                { ByteString("电流 j").FromUtf8(), EMVIEW_J },
                { ByteString("电场+磁场").FromUtf8(), EMVIEW_E_B },
                { ByteString("电场+磁感线").FromUtf8(), EMVIEW_E_B_LINES },
                { ByteString("电场+磁场+电流").FromUtf8(), EMVIEW_E_B_J },
                { ByteString("电场+磁感线+电流").FromUtf8(), EMVIEW_E_B_LINES_J },
                { ByteString("磁场强度 H").FromUtf8(), EMVIEW_H },
                { ByteString("磁化强度 M").FromUtf8(), EMVIEW_M },
                { ByteString("材料类型").FromUtf8(), EMVIEW_TYPE },
                { ByteString("矢量位 A").FromUtf8(), EMVIEW_A },
                { ByteString("坡印廷矢量").FromUtf8(), EMVIEW_POYNTING },
                { ByteString("能量密度").FromUtf8(), EMVIEW_ENERGY },
                { ByteString("坡印廷+能量").FromUtf8(), EMVIEW_POYNTING_ENERGY },
                { ByteString("受力").FromUtf8(), EMVIEW_FORCE },
                { ByteString("等效电流").FromUtf8(), EMVIEW_EFF_CUR },
                { ByteString("磁荷").FromUtf8(), EMVIEW_MAG_CHARGE },
                { ByteString("电场旋度").FromUtf8(), EMVIEW_CURL_E },
                { ByteString("Bx 分量").FromUtf8(), EMVIEW_BX },
                { ByteString("By 分量").FromUtf8(), EMVIEW_BY },
                { ByteString("Hx 分量").FromUtf8(), EMVIEW_HX },
                { ByteString("Hy 分量").FromUtf8(), EMVIEW_HY },
        }, [this] {
                c->SetEMViewMode(emViewMode->GetOption().second);
        });
        emSourceMode = addDropDown(ByteString("电磁场源模式").FromUtf8(), {
                { ByteString("无源").FromUtf8(), EMSRC_NONE },
                { ByteString("单源单频").FromUtf8(), EMSRC_1S1F },
                { ByteString("单源双频").FromUtf8(), EMSRC_1S2F },
                { ByteString("双源单频").FromUtf8(), EMSRC_2S1F },
                { ByteString("双源双频").FromUtf8(), EMSRC_2S2F },
                { ByteString("三源单频").FromUtf8(), EMSRC_3S1F },
                { ByteString("四源单频").FromUtf8(), EMSRC_4S1F },
                { ByteString("单源波包").FromUtf8(), EMSRC_1S1F_PACKET },
                { ByteString("平面波单频").FromUtf8(), EMSRC_1S1F_PLANE },
                { ByteString("平面波双频").FromUtf8(), EMSRC_1S2F_PLANE },
                { ByteString("双平面波").FromUtf8(), EMSRC_2S1F_PLANE },
                { ByteString("双平面双频").FromUtf8(), EMSRC_2S2F_PLANE },
                { ByteString("平面波包").FromUtf8(), EMSRC_1S1F_PLANE_PACKET },
        }, [this] {
                c->SetEMSourceMode(emSourceMode->GetOption().second);
        });
        emFrequency = addTextboxWithPreview(ByteString("电磁波频率 (1~40)").FromUtf8(), false, [this](String value, bool defocus) {
                UpdateEMFrequency(value, defocus);
        }).first;
        emBrightness = addTextboxWithPreview(ByteString("电磁场显示亮度 (1~500)").FromUtf8(), false, [this](String value, bool defocus) {
                UpdateEMBrightness(value, defocus);
        }).first;
        emAux = addTextboxWithPreview(ByteString("辅助参数: 相位/第二频率 (1~40)").FromUtf8(), false, [this](String value, bool defocus) {
                UpdateEMAux(value, defocus);
        }).first;
        emSpeed = addDropDown(ByteString("电磁场模拟速度").FromUtf8(), {
                { "0.5x", 0 },
                { "1x", 1 },
                { "2x", 2 },
        }, [this] {
                c->SetEMSpeed(emSpeed->GetOption().second);
        });
        {
                auto *emClearButton = new ui::Button(ui::Point(10, currentY), ui::Point(90, 16), ByteString("清除电磁场").FromUtf8());
                emClearButton->SetActionCallback({ [this] {
                        c->ClearEMField();
                } });
                scrollPanel->AddChild(emClearButton);
                auto *emClearOvButton = new ui::Button(ui::Point(Size.X - 108, currentY), ui::Point(90, 16), ByteString("清除EM覆盖").FromUtf8());
                emClearOvButton->SetActionCallback({ [this] {
                        c->ClearEMOverrides();
                } });
                scrollPanel->AddChild(emClearOvButton);
                currentY += 26;
        }

        currentY += 4;
        if constexpr (ALLOW_DATA_FOLDER)
        {
                auto *dataFolderButton = new ui::Button(ui::Point(10, currentY), ui::Point(90, 16), ByteString("打开数据目录").FromUtf8());
                dataFolderButton->SetActionCallback({ [] {
                        ByteString cwd = Platform::GetCwd();
                        if (!cwd.empty())
                        {
                                Platform::OpenURI(cwd);
                        }
                        else
                        {
                                std::cerr << "无法打开数据目录: Platform::GetCwd(...) 失败" << std::endl;
                        }
                } });
                scrollPanel->AddChild(dataFolderButton);
                if constexpr (SHARED_DATA_FOLDER)
                {
                        auto *migrationButton = new ui::Button(ui::Point(Size.X - 178, currentY), ui::Point(163, 16), ByteString("迁移至用户共享数据目录").FromUtf8());
                        migrationButton->SetActionCallback({ [] {
                                ByteString from = Platform::originalCwd;
                                ByteString to = Platform::sharedCwd;
                                new ConfirmPrompt(ByteString("进行迁移?").FromUtf8(), ByteString("这将从中迁移所有图章,沙盘和脚本\n\bt").FromUtf8() + from.FromUtf8() + ByteString("\bw\n到共享数据目录\n\bt").FromUtf8() + to.FromUtf8() + "\bw\n\n" + ByteString("已经存在的文件不会被覆盖").FromUtf8(), { [from, to]() {
                                        String ret = Client::Ref().DoMigration(from, to);
                                        new InformationMessage(ByteString("迁移完成").FromUtf8(), ret, false);
                                } });
                        } });
                        scrollPanel->AddChild(migrationButton);
                }
                currentY += 26;
        }
        String autoStartupRequestNote = ByteString("启动时完成").FromUtf8();
        if (!IGNORE_UPDATES)
        {
                autoStartupRequestNote += ByteString(", 保持检查更新").FromUtf8();
        }
        autoStartupRequest = addCheckbox(0, ByteString("获取当天的消息和通知").FromUtf8(), autoStartupRequestNote, [this] {
                auto checked = autoStartupRequest->GetChecked();
                if (checked)
                {
                        Client::Ref().BeginStartupRequest();
                }
                c->SetAutoStartupRequest(checked);
        });
        auto *doStartupRequest = new ui::Button(ui::Point(10, currentY), ui::Point(90, 16), ByteString("立刻获取").FromUtf8());
        doStartupRequest->SetActionCallback({ [] {
                Client::Ref().BeginStartupRequest();
        } });
        scrollPanel->AddChild(doStartupRequest);
        startupRequestStatus = addLabel(5, "");
        UpdateStartupRequestStatus();
        currentY += 13;
        redirectStd = addCheckbox(0, ByteString("将错误等信息保存至文件").FromUtf8(), ByteString("便于开发者排查问题").FromUtf8(), [this] {
                c->SetRedirectStd(redirectStd->GetChecked());
        });

        {
                addSeparator();

                auto *creditsButton = new ui::Button(ui::Point(10, currentY), ui::Point(90, 16), "贡献者");
                creditsButton->SetActionCallback({ [] {
                        auto *credits = new Credits();
                        ui::Engine::Ref().ShowWindow(credits);
                } });
                scrollPanel->AddChild(creditsButton);

                addLabel(5,  ByteString(" - 查看谁为TPT做出了贡献").FromUtf8());
                currentY += 13;
        }


        {
                ui::Button *ok = new ui::Button(ui::Point(0, Size.Y-16), ui::Point(Size.X, 16), "OK");
                ok->SetActionCallback({ [this] {
                        c->Exit();
                } });
                AddComponent(ok);
                SetCancelButton(ok);
                SetOkayButton(ok);
        }
        scrollPanel->InnerSize = ui::Point(Size.X, currentY);
}

void OptionsView::UpdateAmbientAirTempPreview(float airTemp, bool isValid)
{
        if (isValid)
        {
                ambientAirTempPreview->Appearance.BackgroundInactive = HeatToColour(airTemp, MIN_TEMP, MAX_TEMP).WithAlpha(0xFF);
                ambientAirTempPreview->SetText("");
        }
        else
        {
                ambientAirTempPreview->Appearance.BackgroundInactive = ui::Colour(0, 0, 0);
                ambientAirTempPreview->SetText("?");
        }
        ambientAirTempPreview->Appearance.BackgroundHover = ambientAirTempPreview->Appearance.BackgroundInactive;
}

void OptionsView::AmbientAirTempToTextBox(float airTemp)
{
        StringBuilder sb;
        sb << Format::Precision(2);
        format::RenderTemperature(sb, airTemp, TempScale(temperatureScale->GetOption().second));
        ambientAirTemp->SetText(sb.Build());
}

void OptionsView::UpdateEdgePressurePreview(float edgePres, bool isValid)
{
        if (isValid)
        {
                edgePressurePreview->Appearance.BackgroundInactive = PressureToColour(edgePres).WithAlpha(0xFF);
                edgePressurePreview->SetText("");
        }
        else
        {
                edgePressurePreview->Appearance.BackgroundInactive = ui::Colour(0, 0, 0);
                edgePressurePreview->SetText("?");
        }
        edgePressurePreview->Appearance.BackgroundHover = edgePressurePreview->Appearance.BackgroundInactive;
}

void OptionsView::EdgePressureToTextBox(float edgePres)
{
        StringBuilder sb;
        sb << Format::Precision(2) << edgePres;
        edgePressure->SetText(sb.Build());
}

void OptionsView::VorticityCoeffToTextBox(float vorticity)
{
        StringBuilder sb;
        sb << Format::Precision(2) << vorticity;
        vorticityCoeff->SetText(sb.Build());
}

void OptionsView::UpdateStartupRequestStatus()
{
        switch (Client::Ref().GetStartupRequestStatus())
        {
        case Client::StartupRequestStatus::notYetDone:
                startupRequestStatus->SetText("\bg - 尚未获取");
                break;

        case Client::StartupRequestStatus::inProgress:
                startupRequestStatus->SetText("\bg - 进行中...");
                break;

        case Client::StartupRequestStatus::succeeded:
                startupRequestStatus->SetText(String::Build("\bg - 完成, ", Client::Ref().GetServerNotifications().size(), " 个通知已获取"));
                break;

        case Client::StartupRequestStatus::failed:
                {
                        auto error = Client::Ref().GetStartupRequestError();
                        if (!error)
                        {
                                error = "???";
                        }
                        startupRequestStatus->SetText("\bg - 失败: " + error->FromUtf8());
                }
                break;
        }
}

static void UpdateSettingFromString(String pres, bool isDefocus, float minVale, float maxValue, float defaultValue, auto parse, auto toTextbox, auto apply)
{
        // Parse value and determine validity
        float value = 0;
        bool isValid;
        try
        {
                value = parse(pres);
                isValid = true;
        }
        catch (const std::exception &ex)
        {
                isValid = false;
        }

        // While defocusing, correct out of range values and empty textboxes
        if (isDefocus)
        {
                if (pres.empty())
                {
                        isValid = true;
                        value = defaultValue;
                }
                else if (!isValid)
                        return;
                else if (value < minVale)
                        value = minVale;
                else if (value > maxValue)
                        value = maxValue;

                toTextbox(value);
        }
        // Out of range values are invalid, preview should go away
        else if (isValid && (value < minVale || value > maxValue))
                isValid = false;

        // If valid, apply
        apply(value, isValid);
}

void OptionsView::UpdateAirTemp(String temp, bool isDefocus)
{
        UpdateSettingFromString(temp, isDefocus, MIN_TEMP, MAX_TEMP, float(R_TEMP) + 273.15f, [this](const String &temp) {
                return format::StringToTemperature(temp, TempScale(temperatureScale->GetOption().second));
        }, [this](float airTemp) {
                AmbientAirTempToTextBox(airTemp);
        }, [this](float airTemp, bool isValid) {
                if (isValid)
                {
                        c->SetAmbientAirTemperature(airTemp);
                }
                UpdateAmbientAirTempPreview(airTemp, isValid);
        });
}

void OptionsView::UpdateEdgePressure(String pres, bool isDefocus)
{
        UpdateSettingFromString(pres, isDefocus, MIN_PRESSURE, MAX_PRESSURE, 0.f, [](const String &pres) {
                return pres.ToNumber<float>();
        }, [this](float edgePres) {
                EdgePressureToTextBox(edgePres);
        }, [this](float edgePres, bool isValid) {
                if (isValid)
                {
                        c->SetEdgePressure(edgePres);
                }
                UpdateEdgePressurePreview(edgePres, isValid);
        });
}

void OptionsView::UpdateVorticityCoeff(String vort, bool isDefocus)
{
        UpdateSettingFromString(vort, isDefocus, 0.f, 1.f, 0.1f, [](const String &vort) {
                return vort.ToNumber<float>();
        }, [this](float vorticity) {
                VorticityCoeffToTextBox(vorticity);
        }, [this](float vorticity, bool isValid) {
                if (isValid)
                {
                        c->SetVorticityCoeff(vorticity);
                }
        });
}

// TM-mode EM field settings (EMWave2 port)

void OptionsView::EMFrequencyToTextBox(float value)
{
        StringBuilder sb;
        sb << Format::Precision(2) << value;
        emFrequency->SetText(sb.Build());
}

void OptionsView::EMBrightnessToTextBox(int value)
{
        StringBuilder sb;
        sb << value;
        emBrightness->SetText(sb.Build());
}

void OptionsView::EMAuxToTextBox(float value)
{
        StringBuilder sb;
        sb << Format::Precision(2) << value;
        emAux->SetText(sb.Build());
}

void OptionsView::UpdateEMFrequency(String value, bool isDefocus)
{
        UpdateSettingFromString(value, isDefocus, 1.f, 40.f, 10.f, [](const String &pres) {
                return pres.ToNumber<float>();
        }, [this](float frequency) {
                EMFrequencyToTextBox(frequency);
        }, [this](float frequency, bool isValid) {
                if (isValid)
                {
                        c->SetEMFrequency(frequency);
                }
        });
}

void OptionsView::UpdateEMBrightness(String value, bool isDefocus)
{
        UpdateSettingFromString(value, isDefocus, 1.f, 500.f, 100.f, [](const String &pres) {
                return pres.ToNumber<float>();
        }, [this](float brightness) {
                EMBrightnessToTextBox(int(brightness));
        }, [this](float brightness, bool isValid) {
                if (isValid)
                {
                        c->SetEMBrightness(int(brightness));
                }
        });
}

void OptionsView::UpdateEMAux(String value, bool isDefocus)
{
        UpdateSettingFromString(value, isDefocus, 1.f, 40.f, 1.f, [](const String &pres) {
                return pres.ToNumber<float>();
        }, [this](float aux) {
                EMAuxToTextBox(aux);
        }, [this](float aux, bool isValid) {
                if (isValid)
                {
                        c->SetEMAux(aux);
                }
        });
}

void OptionsView::NotifySettingsChanged(OptionsModel * sender)
{
        temperatureScale->SetOption(sender->GetTemperatureScale()); // has to happen before AmbientAirTempToTextBox is called
        heatSimulation->SetChecked(sender->GetHeatSimulation());
        ambientHeatSimulation->SetChecked(sender->GetAmbientHeatSimulation());
        // BUGFIX: language is only created when FORCE_WINDOW_FRAME_OPS != forceWindowFrameOpsHandheld
        // (e.g. it is missing on handheld/android builds); dereferencing it unconditionally here
        // caused a "Memory read/write error" crash whenever the settings window was opened there.
        if (language)
        {
                language->SetOption(sender->GetLanguage());
        }
        newtonianGravity->SetChecked(sender->GetNewtonianGravity());
        waterEqualisation->SetChecked(sender->GetWaterEqualisation());
        airMode->SetOption(sender->GetAirMode());
        // Initialize air temp and preview only when the options menu is opened, and not when user is actively editing the textbox
        if (!ambientAirTemp->IsFocused())
        {
                float airTemp = sender->GetAmbientAirTemperature();
                UpdateAmbientAirTempPreview(airTemp, true);
                AmbientAirTempToTextBox(airTemp);
        }
        if (!edgePressure->IsFocused())
        {
                float pres = sender->GetEdgePressure();
                UpdateEdgePressurePreview(pres, true);
                EdgePressureToTextBox(pres);
        }
        // Same for vorticity
        if (!vorticityCoeff->IsFocused())
        {
                VorticityCoeffToTextBox(sender->GetVorticityCoeff());
        }
        convectionMode->SetOption(sender->GetConvectionMode());
        edgeVelocityX = sender->GetEdgeVelocityX();
        edgeVelocityY = sender->GetEdgeVelocityY();
        gravityMode->SetOption(sender->GetGravityMode());
        customGravityX = sender->GetCustomGravityX();
        customGravityY = sender->GetCustomGravityY();
        decoSpace->SetOption(sender->GetDecoSpace());
        edgeMode->SetOption(sender->GetEdgeMode());
        if (scale)
        {
                scale->SetOption(sender->GetScale());
        }
        if (resizable)
        {
                resizable->SetChecked(sender->GetResizable());
        }
        if (fullscreen)
        {
                fullscreen->SetChecked(sender->GetFullscreen());
        }
        if (changeResolution)
        {
                changeResolution->SetChecked(sender->GetChangeResolution());
        }
        if (forceIntegerScaling)
        {
                forceIntegerScaling->SetChecked(sender->GetForceIntegerScaling());
        }
        if (blurryScaling)
        {
                blurryScaling->SetChecked(sender->GetBlurryScaling());
        }
        if (fastquit)
        {
                fastquit->SetChecked(sender->GetFastQuit());
        }
        if (globalQuit)
        {
                globalQuit->SetChecked(sender->GetGlobalQuit());
        }
        if (nativeClipoard)
        {
                nativeClipoard->SetChecked(sender->GetNativeClipoard());
        }
        showAvatars->SetChecked(sender->GetShowAvatars());
        mouseClickRequired->SetChecked(sender->GetMouseClickRequired());
        includePressure->SetChecked(sender->GetIncludePressure());
        perfectCircle->SetChecked(sender->GetPerfectCircle());
        graveExitsConsole->SetChecked(sender->GetGraveExitsConsole());
        threadedRendering->SetChecked(sender->GetThreadedRendering());
        momentumScroll->SetChecked(sender->GetMomentumScroll());
        redirectStd->SetChecked(sender->GetRedirectStd());
        autoStartupRequest->SetChecked(sender->GetAutoStartupRequest());
        // TM-mode EM field settings (EMWave2 port); each component is guarded,
        // exactly like the language dropdown fix above, so a missing component
        // can never crash the settings window
        if (emEnabled)
        {
                emEnabled->SetChecked(sender->GetEMEnabled());
        }
        if (emCellSize)
        {
                emCellSize->SetOption(sender->GetEMCellSize());
        }
        if (emSourceMode)
        {
                emSourceMode->SetOption(sender->GetEMSourceMode());
        }
        if (emViewMode)
        {
                emViewMode->SetOption(sender->GetEMViewMode());
        }
        if (emSpeed)
        {
                emSpeed->SetOption(sender->GetEMSpeed());
        }
        if (emFrequency && !emFrequency->IsFocused())
        {
                EMFrequencyToTextBox(sender->GetEMFrequency());
        }
        if (emBrightness && !emBrightness->IsFocused())
        {
                EMBrightnessToTextBox(sender->GetEMBrightness());
        }
        if (emAux && !emAux->IsFocused())
        {
                EMAuxToTextBox(sender->GetEMAux());
        }
}

void OptionsView::AttachController(OptionsController * c_)
{
        c = c_;
}

void OptionsView::OnTick()
{
#if EMFIELD_DEBUG
        // debug block: the EM self test opens this window and can't interact with
        // it headlessly, so close it again after a short delay to let the test
        // continue; compiled out of normal builds
        {
                static int emDebugTicks = 0;
                if (++emDebugTicks >= 30 && c)
                {
                        std::cout << "[EMSELFTEST] options window ticked " << emDebugTicks
                                  << " frames, closing it again" << std::endl;
                        emDebugTicks = 0;
                        c->Exit();
                }
        }
#endif
        UpdateStartupRequestStatus();
}

void OptionsView::OnDraw()
{
        Graphics * g = GetGraphics();
        g->DrawFilledRect(RectSized(Position - Vec2{ 1, 1 }, Size + Vec2{ 2, 2 }), 0x000000_rgb);
        g->DrawRect(RectSized(Position, Size), 0xFFFFFF_rgb);
}

void OptionsView::OnTryExit(ExitMethod method)
{
        c->Exit();
}
