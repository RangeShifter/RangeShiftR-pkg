object frmSpecies: TfrmSpecies
  Left = 0
  Top = 0
  Width = 1200
  Height = 641
  AutoScroll = True
  Caption = 'Species Parameters'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poOwnerFormCenter
  PixelsPerInch = 96
  TextHeight = 13
  object PageControl1: TPageControl
    Left = 0
    Top = 0
    Width = 1184
    Height = 569
    ActivePage = TSPop
    Align = alTop
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 0
    object TSPop: TTabSheet
      Caption = 'Population dynamics'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
      object PanelK: TPanel
        Left = 633
        Top = 11
        Width = 360
        Height = 253
        Alignment = taLeftJustify
        BevelInner = bvLowered
        TabOrder = 9
        VerticalAlignment = taAlignTop
        object SGhabLabel: TLabel
          Left = 14
          Top = 50
          Width = 260
          Height = 16
          Caption = 'Habitat-dependent carrying capacities: '
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentFont = False
          WordWrap = True
        end
        object edtKLabel: TLabel
          Left = 84
          Top = 23
          Width = 200
          Height = 16
          Caption = 'inds / ha (assuming 100% quality) '
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
        end
        object SGhab: TStringGrid
          Left = 14
          Top = 72
          Width = 211
          Height = 169
          BiDiMode = bdLeftToRight
          ColCount = 2
          DefaultColWidth = 90
          DoubleBuffered = False
          RowCount = 12
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 1
          ColWidths = (
            90
            90)
          RowHeights = (
            24
            24
            24
            24
            24
            24
            24
            24
            24
            24
            24
            24)
        end
        object edtK: TLabeledEdit
          Left = 34
          Top = 20
          Width = 44
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 16
          EditLabel.Height = 16
          EditLabel.Caption = 'K: '
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpLeft
          LabelSpacing = 2
          ParentFont = False
          TabOrder = 0
          Text = '10.0'
        end
      end
      object RGReproduction: TRadioGroup
        Left = 3
        Top = 61
        Width = 182
        Height = 100
        Ctl3D = True
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ItemIndex = 0
        Items.Strings = (
          'Asexual / Only females model'
          'Simple sexual model'
          'Complex sexual model')
        ParentCtl3D = False
        ParentFont = False
        TabOrder = 1
        WordWrap = True
        OnClick = RGReproductionClick
      end
      object edtR: TLabeledEdit
        Left = 344
        Top = 83
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 130
        EditLabel.Height = 32
        EditLabel.Caption = 'Rmax              (intrinsic growth rate):'
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        EditLabel.WordWrap = True
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        LabelSpacing = 2
        ParentFont = False
        TabOrder = 3
        Text = '1.5'
      end
      object CBStageModel: TCheckBox
        Left = 3
        Top = 3
        Width = 222
        Height = 55
        Caption = 'Overlapping generations / Stage-structured model '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        WordWrap = True
        OnClick = CBStageModelClick
      end
      object PanelStage: TPanel
        Left = 4
        Top = 270
        Width = 989
        Height = 265
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Caption = 'Stage-structured population model '
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 10
        VerticalAlignment = taAlignTop
        Visible = False
        object transMatrixLabel: TLabel
          Left = 9
          Top = 59
          Width = 117
          Height = 16
          Caption = 'Transition Matrix  '
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentFont = False
        end
        object Bevel1: TBevel
          Left = 511
          Top = 32
          Width = 13
          Height = 225
          Shape = bsRightLine
        end
        object MinAgesLabel: TLabel
          Left = 393
          Top = 59
          Width = 96
          Height = 16
          Caption = 'Minimum Ages '
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentFont = False
        end
        object Bevel2: TBevel
          Left = 681
          Top = 32
          Width = 13
          Height = 225
          Shape = bsRightLine
        end
        object SurvivalLabel: TLabel
          Left = 530
          Top = 31
          Width = 145
          Height = 16
          Caption = 'Scheduling of survival '
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentFont = False
        end
        object edtNstages: TLabeledEdit
          Left = 97
          Top = 24
          Width = 44
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 82
          EditLabel.Height = 16
          EditLabel.Caption = 'Nr. of stages: '
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpLeft
          LabelSpacing = 2
          ParentFont = False
          TabOrder = 0
          Text = '2'
          OnExit = edtNstagesExit
        end
        object transMatrix: TStringGrid
          Left = 9
          Top = 81
          Width = 378
          Height = 175
          BiDiMode = bdLeftToRight
          ColCount = 3
          DefaultColWidth = 50
          DoubleBuffered = False
          RowCount = 3
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 2
          ColWidths = (
            50
            50
            50)
          RowHeights = (
            24
            24
            24)
        end
        object edtMaxAge: TLabeledEdit
          Left = 225
          Top = 24
          Width = 44
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 61
          EditLabel.Height = 16
          EditLabel.Caption = 'Max. age: '
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpLeft
          LabelSpacing = 2
          ParentFont = False
          TabOrder = 1
          Text = '1000'
        end
        object MinAges: TStringGrid
          Left = 393
          Top = 81
          Width = 124
          Height = 175
          BiDiMode = bdLeftToRight
          ColCount = 2
          DefaultColWidth = 50
          DoubleBuffered = False
          RowCount = 2
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 3
          OnExit = MinAgesExit
          ColWidths = (
            50
            50)
          RowHeights = (
            24
            24)
        end
        object GroupBox1: TGroupBox
          Left = 700
          Top = 32
          Width = 147
          Height = 105
          Caption = 'Density dependence '
          TabOrder = 5
          object CBFecundity: TCheckBox
            Left = 16
            Top = 24
            Width = 97
            Height = 17
            Caption = 'Fecundity '
            Checked = True
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            State = cbChecked
            TabOrder = 0
            OnClick = CBFecundityClick
          end
          object CBDevelopment: TCheckBox
            Left = 16
            Top = 47
            Width = 97
            Height = 17
            Caption = 'Development '
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 1
            OnClick = CBDevelopmentClick
          end
          object CBSurvival: TCheckBox
            Left = 16
            Top = 70
            Width = 97
            Height = 17
            Caption = 'Survival '
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 2
            OnClick = CBSurvivalClick
          end
        end
        object GroupBox2: TGroupBox
          Left = 700
          Top = 143
          Width = 147
          Height = 105
          Caption = 'Stages'#39' weights '
          TabOrder = 8
          object CBweightFec: TCheckBox
            Left = 16
            Top = 24
            Width = 97
            Height = 17
            Caption = 'Fecundity '
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 0
            OnClick = CBweightFecClick
          end
          object CBweightDev: TCheckBox
            Left = 16
            Top = 47
            Width = 97
            Height = 17
            Caption = 'Development '
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 1
            OnClick = CBweightDevClick
          end
          object CBweightSurv: TCheckBox
            Left = 16
            Top = 70
            Width = 97
            Height = 17
            Caption = 'Survival '
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 2
            OnClick = CBweightSurvClick
          end
        end
        object BtnSetCoeff: TButton
          Left = 715
          Top = 236
          Width = 142
          Height = 26
          Caption = 'Set weights'
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 9
          WordWrap = True
          OnClick = BtnSetCoeffClick
        end
        object edtDevCoeff: TLabeledEdit
          Left = 853
          Top = 72
          Width = 44
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 78
          EditLabel.Height = 32
          EditLabel.Caption = 'Development coeff.'
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          EditLabel.WordWrap = True
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpRight
          LabelSpacing = 2
          ParentFont = False
          TabOrder = 6
          Text = '1.0'
          Visible = False
        end
        object edtSurvCoeff: TLabeledEdit
          Left = 853
          Top = 102
          Width = 44
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 81
          EditLabel.Height = 16
          EditLabel.Caption = 'Survival coeff.'
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpRight
          LabelSpacing = 2
          ParentFont = False
          TabOrder = 7
          Text = '1.0'
          Visible = False
        end
        object RGSurvival: TRadioGroup
          Left = 530
          Top = 53
          Width = 160
          Height = 116
          Ctl3D = True
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ItemIndex = 1
          Items.Strings = (
            'At reproduction'
            'Between reproductive events'
            'Annually')
          ParentCtl3D = False
          ParentFont = False
          TabOrder = 4
          WordWrap = True
        end
      end
      object edtSexRatio: TLabeledEdit
        Left = 344
        Top = 137
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 121
        EditLabel.Height = 16
        EditLabel.Caption = 'Proportion of males: '
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        LabelSpacing = 2
        ParentFont = False
        TabOrder = 5
        Text = '0.5'
        Visible = False
      end
      object edtHarem: TLabeledEdit
        Left = 568
        Top = 137
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 153
        EditLabel.Height = 16
        EditLabel.Caption = 'h (maximum harem size): '
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        LabelSpacing = 2
        ParentFont = False
        TabOrder = 6
        Text = '1'
        Visible = False
      end
      object edtC: TLabeledEdit
        Left = 568
        Top = 83
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 164
        EditLabel.Height = 16
        EditLabel.Caption = 'bc (competition coefficient): '
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        LabelSpacing = 2
        ParentFont = False
        TabOrder = 4
        Text = '1.0'
      end
      object edtPRep: TLabeledEdit
        Left = 568
        Top = 190
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 155
        EditLabel.Height = 16
        EditLabel.Caption = 'Probability of reproducing: '
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        LabelSpacing = 2
        ParentFont = False
        TabOrder = 7
        Text = '1.0'
        Visible = False
      end
      object edtRepInterval: TLabeledEdit
        Left = 568
        Top = 228
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 202
        EditLabel.Height = 32
        EditLabel.Caption = 'Nr. of reproductive seasons before subsequent reproduction: '
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        EditLabel.WordWrap = True
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        LabelSpacing = 2
        ParentFont = False
        TabOrder = 8
        Text = '0'
        Visible = False
      end
      object edtGen: TLabeledEdit
        Left = 417
        Top = 18
        Width = 44
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 189
        EditLabel.Height = 16
        EditLabel.Caption = 'Nr. reproductive seasons / year: '
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentFont = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        LabelPosition = lpLeft
        ParentFont = False
        TabOrder = 2
        Text = '1'
      end
      object RGaction: TRadioGroup
        Left = 3
        Top = 176
        Width = 238
        Height = 76
        Caption = 'Action after population destruction: '
        Ctl3D = True
        DoubleBuffered = False
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ItemIndex = 0
        Items.Strings = (
          'All individuals die'
          'All individuals disperse')
        ParentCtl3D = False
        ParentDoubleBuffered = False
        ParentFont = False
        TabOrder = 11
        Visible = False
      end
    end
    object TSDispersal: TTabSheet
      Caption = 'Dispersal'
      ImageIndex = 1
      OnExit = TSDispersalExit
      object PanelSettProcess: TPanel
        Left = 703
        Top = 3
        Width = 458
        Height = 532
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 2
        VerticalAlignment = taAlignTop
        Visible = False
        object SettProcessLabel: TLabel
          Left = 8
          Top = 2
          Width = 240
          Height = 16
          Caption = 'SETTLEMENT - Movement Processes   '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object PanelMateDD: TPanel
          Left = 5
          Top = 80
          Width = 444
          Height = 129
          BevelOuter = bvNone
          TabOrder = 5
          object CBFindMate: TCheckBox
            Left = 10
            Top = 20
            Width = 222
            Height = 17
            Caption = 'Find mate'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 1
          end
          object CBDensDepSettle: TCheckBox
            Left = 10
            Top = 50
            Width = 222
            Height = 17
            Caption = 'Density dependence'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 0
            OnClick = CBDensDepSettleClick
          end
          object CBIndVarSettle: TCheckBox
            Left = 10
            Top = 80
            Width = 222
            Height = 17
            Caption = 'Individual variability'
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 2
            OnClick = CBIndVarSettleClick
          end
        end
        object CBSexSettMovt: TCheckBox
          Left = 280
          Top = 25
          Width = 157
          Height = 17
          Caption = 'Sex dependent'
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold]
          ParentFont = False
          TabOrder = 6
          Visible = False
          OnClick = CBSexSettMovtClick
        end
        object CBStageSettMovt: TCheckBox
          Left = 280
          Top = 48
          Width = 141
          Height = 17
          Caption = 'Stage dependent'
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold]
          ParentFont = False
          TabOrder = 7
          Visible = False
          OnClick = CBStageSettMovtClick
        end
        object edtMinSteps: TLabeledEdit
          Left = 111
          Top = 48
          Width = 40
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 97
          EditLabel.Height = 16
          EditLabel.Caption = 'Min. nr. of steps '
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = [fsItalic]
          EditLabel.ParentFont = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpLeft
          ParentFont = False
          TabOrder = 8
          Text = '0'
        end
        object PanelSettleDD: TPanel
          Left = 8
          Top = 215
          Width = 441
          Height = 150
          Alignment = taLeftJustify
          BevelOuter = bvNone
          Caption = 'Density dependence parameters:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 0
          VerticalAlignment = taAlignTop
          Visible = False
          DesignSize = (
            441
            150)
          object edtSettS0: TLabeledEdit
            Left = 110
            Top = 50
            Width = 40
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 24
            EditLabel.Height = 16
            EditLabel.Caption = 'S0: '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 0
            Text = '1'
          end
          object edtSettAlpha: TLabeledEdit
            Left = 110
            Top = 80
            Width = 40
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 49
            EditLabel.Height = 16
            EditLabel.Caption = 'AlphaS: '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            EditLabel.ShowAccelChar = False
            EditLabel.WordWrap = True
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 4
            Text = '0'
          end
          object edtSettBeta: TLabeledEdit
            Left = 110
            Top = 110
            Width = 40
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 42
            EditLabel.Height = 16
            EditLabel.Caption = 'BetaS: '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 5
            Text = '1'
          end
          object edtSettS0Mean: TLabeledEdit
            Left = 190
            Top = 50
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 36
            EditLabel.Height = 16
            EditLabel.Caption = 'Mean:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 1
            Text = '0.0'
            Visible = False
          end
          object edtSettS0SD: TLabeledEdit
            Left = 270
            Top = 50
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 24
            EditLabel.Height = 16
            EditLabel.Caption = 'S.d:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 2
            Text = '0.0'
            Visible = False
          end
          object edtSettAlphaMean: TEdit
            Left = 190
            Top = 80
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 6
            Text = '0.0'
            Visible = False
          end
          object edtSettBetaMean: TEdit
            Left = 190
            Top = 110
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 9
            Text = '0.0'
            Visible = False
          end
          object edtSettAlphaSD: TEdit
            Left = 270
            Top = 80
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 8
            Text = '0.0'
            Visible = False
          end
          object edtSettBetaSD: TEdit
            Left = 270
            Top = 110
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 10
            Text = '0.0'
            Visible = False
          end
          object edtSettS0Scale: TLabeledEdit
            Left = 350
            Top = 50
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.AlignWithMargins = True
            EditLabel.Width = 45
            EditLabel.Height = 32
            EditLabel.BiDiMode = bdLeftToRight
            EditLabel.Caption = 'Scaling factor:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentBiDiMode = False
            EditLabel.ParentFont = False
            EditLabel.WordWrap = True
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 3
            Text = '0.0'
            Visible = False
          end
          object edtSettAlphaScale: TEdit
            Left = 350
            Top = 80
            Width = 50
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 7
            Text = '0.0'
            Visible = False
          end
          object edtSettBetaScale: TEdit
            Left = 350
            Top = 110
            Width = 50
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 11
            Text = '0.0'
            Visible = False
          end
        end
        object RGSteps: TRadioGroup
          Left = 13
          Top = 375
          Width = 268
          Height = 68
          Caption = 'If not settled, move until...'
          DoubleBuffered = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ItemIndex = 1
          Items.Strings = (
            'Maximum nr. of steps '
            'Only per-step mortality')
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 3
          WordWrap = True
          OnClick = RGStepsClick
        end
        object edtNsteps: TLabeledEdit
          Left = 176
          Top = 448
          Width = 40
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 105
          EditLabel.Height = 16
          EditLabel.Caption = 'Max. nr. of steps: '
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpLeft
          ParentFont = False
          TabOrder = 1
          Text = '0'
          Visible = False
        end
        object edtMaxStepYear: TLabeledEdit
          Left = 176
          Top = 480
          Width = 40
          Height = 24
          Alignment = taRightJustify
          EditLabel.Width = 157
          EditLabel.Height = 16
          EditLabel.Caption = 'Max. nr. of steps per year: '
          EditLabel.Font.Charset = DEFAULT_CHARSET
          EditLabel.Font.Color = clWindowText
          EditLabel.Font.Height = -13
          EditLabel.Font.Name = 'Tahoma'
          EditLabel.Font.Style = []
          EditLabel.ParentFont = False
          EditLabel.WordWrap = True
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          LabelPosition = lpLeft
          ParentFont = False
          TabOrder = 2
          Text = '0'
        end
        object Memo2: TMemo
          Left = 232
          Top = 480
          Width = 217
          Height = 46
          BevelKind = bkFlat
          BorderStyle = bsNone
          Color = clInfoBk
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -11
          Font.Name = 'Tahoma'
          Font.Style = [fsItalic]
          Lines.Strings = (
            'If zero, every individual completes the '
            'dispersal phase in one year (between '
            'two successive reproduction phases)')
          ParentFont = False
          ReadOnly = True
          TabOrder = 4
          Visible = False
        end
      end
      object PanelSettKernels: TPanel
        Left = 703
        Top = 3
        Width = 330
        Height = 245
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 3
        VerticalAlignment = taAlignTop
        object SettKernelsLabel: TLabel
          Left = 8
          Top = 3
          Width = 213
          Height = 16
          Caption = 'SETTLEMENT - Dispersal Kernels   '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object RGSettKern: TRadioGroup
          Left = 8
          Top = 55
          Width = 313
          Height = 154
          Caption = 'If the arrival cell is unsuitable:'
          Ctl3D = True
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ItemIndex = 0
          Items.Strings = (
            'Die'
            'Wait'
            'SET IN CODE'
            'SET IN CODE')
          ParentCtl3D = False
          ParentFont = False
          TabOrder = 2
          WordWrap = True
          OnClick = RGSettKernClick
        end
        object CBSexSettKern: TCheckBox
          Left = 12
          Top = 25
          Width = 173
          Height = 17
          Caption = 'Sex dependent'
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold]
          ParentFont = False
          TabOrder = 0
          Visible = False
          OnClick = CBSexSettKernClick
        end
        object CBStageSettKern: TCheckBox
          Left = 148
          Top = 25
          Width = 173
          Height = 17
          Caption = 'Stage dependent'
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold]
          ParentFont = False
          TabOrder = 1
          Visible = False
          OnClick = CBStageSettKernClick
        end
        object CBSettKernMate: TCheckBox
          Left = 12
          Top = 215
          Width = 173
          Height = 17
          Caption = '+ mating requirements  '
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 3
          Visible = False
        end
      end
      object PanelEmigration: TPanel
        Left = 3
        Top = 3
        Width = 257
        Height = 532
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 0
        VerticalAlignment = taAlignTop
        DesignSize = (
          257
          532)
        object EmigrationLabel: TLabel
          Left = 8
          Top = 2
          Width = 87
          Height = 16
          Caption = 'EMIGRATION  '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object PanelEP: TPanel
          Left = 0
          Top = 246
          Width = 249
          Height = 100
          BevelEdges = []
          BevelOuter = bvNone
          TabOrder = 1
          DesignSize = (
            249
            100)
          object LabelDensIndept: TLabel
            Left = 9
            Top = 11
            Width = 161
            Height = 16
            Anchors = []
            AutoSize = False
            Caption = 'Density independent (d):'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            WordWrap = True
            ExplicitTop = 8
          end
          object edtEP: TLabeledEdit
            Left = 48
            Top = 56
            Width = 40
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 15
            EditLabel.Height = 16
            EditLabel.Caption = 'd  '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = [fsItalic]
            EditLabel.ParentFont = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 0
            Text = '0.0'
          end
          object edtEPmean: TLabeledEdit
            Left = 94
            Top = 56
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 36
            EditLabel.Height = 16
            EditLabel.Caption = 'Mean:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 1
            Text = '0.0'
            Visible = False
          end
          object edtEPsd: TLabeledEdit
            Left = 141
            Top = 56
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 28
            EditLabel.Height = 16
            EditLabel.Caption = 'S.d.:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 2
            Text = '0.0'
            Visible = False
          end
          object edtEPscale: TLabeledEdit
            Left = 188
            Top = 56
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.AlignWithMargins = True
            EditLabel.Width = 45
            EditLabel.Height = 32
            EditLabel.BiDiMode = bdLeftToRight
            EditLabel.Caption = 'Scaling factor:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentBiDiMode = False
            EditLabel.ParentFont = False
            EditLabel.WordWrap = True
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 3
            Text = '0.0'
            Visible = False
          end
        end
        object RGEmigProb: TRadioGroup
          Left = 8
          Top = 22
          Width = 160
          Height = 67
          Caption = 'Emigration Probability'
          Color = clBtnHighlight
          Ctl3D = True
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ItemIndex = 0
          Items.Strings = (
            'Density-independent'
            'Density-dependent')
          ParentBackground = False
          ParentColor = False
          ParentCtl3D = False
          ParentFont = False
          TabOrder = 0
          OnClick = RGEmigProbClick
        end
        object PanelDensEmig: TPanel
          Left = 8
          Top = 347
          Width = 249
          Height = 176
          Alignment = taLeftJustify
          Anchors = []
          BevelEdges = []
          BevelOuter = bvNone
          TabOrder = 2
          Visible = False
          DesignSize = (
            249
            176)
          object edtD0Mean: TLabeledEdit
            Left = 86
            Top = 57
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 36
            EditLabel.Height = 16
            EditLabel.Caption = 'Mean:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 1
            Text = '0.0'
            Visible = False
          end
          object edtD0SD: TLabeledEdit
            Left = 133
            Top = 57
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 24
            EditLabel.Height = 16
            EditLabel.Caption = 'S.d:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 2
            Text = '0.0'
            Visible = False
          end
          object edtD0: TLabeledEdit
            Left = 40
            Top = 57
            Width = 40
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            EditLabel.Width = 19
            EditLabel.Height = 16
            EditLabel.Caption = 'D0 '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = [fsItalic]
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 0
            Text = '0.0'
          end
          object edtAlpha: TLabeledEdit
            Left = 40
            Top = 91
            Width = 40
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 11
            EditLabel.Height = 16
            EditLabel.Caption = 'a '
            EditLabel.Font.Charset = GREEK_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Symbol'
            EditLabel.Font.Style = [fsItalic]
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = GREEK_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 4
            Text = '0.0'
          end
          object edtBeta: TLabeledEdit
            Left = 40
            Top = 121
            Width = 40
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            EditLabel.Width = 10
            EditLabel.Height = 16
            EditLabel.Caption = 'b '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Symbol'
            EditLabel.Font.Style = [fsItalic]
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 8
            Text = '0.0'
          end
          object edtAlphaMean: TEdit
            Left = 86
            Top = 91
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 5
            Text = '0.0'
            Visible = False
          end
          object edtBetaMean: TEdit
            Left = 86
            Top = 121
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 9
            Text = '0.0'
            Visible = False
          end
          object edtAlphaSD: TEdit
            Left = 133
            Top = 91
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 6
            Text = '0.0'
            Visible = False
          end
          object edtBetaSD: TEdit
            Left = 133
            Top = 121
            Width = 41
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 10
            Text = '0.0'
            Visible = False
          end
          object edtD0Scale: TLabeledEdit
            Left = 180
            Top = 56
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.AlignWithMargins = True
            EditLabel.Width = 45
            EditLabel.Height = 32
            EditLabel.BiDiMode = bdLeftToRight
            EditLabel.Caption = 'Scaling factor:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentBiDiMode = False
            EditLabel.ParentFont = False
            EditLabel.WordWrap = True
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 3
            Text = '0.0'
            Visible = False
          end
          object edtAlphaScale: TEdit
            Left = 180
            Top = 91
            Width = 50
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 7
            Text = '0.0'
            Visible = False
          end
          object edtBetaScale: TEdit
            Left = 180
            Top = 121
            Width = 50
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 11
            Text = '0.0'
            Visible = False
          end
        end
        object EmigCheckboxPanel: TPanel
          Left = 8
          Top = 95
          Width = 229
          Height = 130
          BevelEdges = []
          BevelOuter = bvNone
          ParentColor = True
          TabOrder = 3
          OnExit = emigCheckboxPanelExit
          object CBFullKernel: TCheckBox
            Left = 8
            Top = 12
            Width = 173
            Height = 17
            Caption = 'Use full kernel'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 0
          end
          object CBSexEmig: TCheckBox
            Left = 8
            Top = 35
            Width = 173
            Height = 17
            Caption = 'Sex dependent'
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 1
            Visible = False
            OnClick = CBSexEmigClick
          end
          object CBStageEmig: TCheckBox
            Left = 8
            Top = 58
            Width = 173
            Height = 17
            Caption = 'Stage dependent'
            Color = clBtnFace
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentColor = False
            ParentFont = False
            TabOrder = 2
            Visible = False
            OnClick = CBStageEmigClick
          end
          object CBIndVarEmig: TCheckBox
            Left = 8
            Top = 81
            Width = 157
            Height = 17
            Caption = 'Individual variability'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 3
            OnClick = CBIndVarEmigClick
          end
          object edtEmigStage: TLabeledEdit
            Left = 115
            Top = 105
            Width = 49
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 105
            EditLabel.Height = 16
            EditLabel.Caption = 'Emigration stage: '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            LabelSpacing = 2
            ParentFont = False
            TabOrder = 4
            Text = '0'
            Visible = False
          end
        end
      end
      object PanelTransfer: TPanel
        Left = 263
        Top = 3
        Width = 434
        Height = 532
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 1
        VerticalAlignment = taAlignTop
        object TransferLabel: TLabel
          Left = 9
          Top = 2
          Width = 76
          Height = 16
          Caption = 'TRANSFER   '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object RGMovements: TRadioGroup
          Left = 9
          Top = 24
          Width = 150
          Height = 67
          Caption = 'Movement Model'
          Ctl3D = True
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ItemIndex = 0
          Items.Strings = (
            'Dispersal kernels'
            'Movement processes')
          ParentCtl3D = False
          ParentFont = False
          TabOrder = 0
          OnClick = RGMovementsClick
        end
        object BtnMovements: TBitBtn
          Left = 165
          Top = 48
          Width = 137
          Height = 25
          Caption = 'Set parameters'
          Enabled = False
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Glyph.Data = {
            36030000424D3603000000000000360000002800000010000000100000000100
            18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
            FF00FFFF00FFFF00FF044906055B09066C0C066C0C055E0A044C06FF00FFFF00
            FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF05600905600908911309B01809
            B31A09B31909B11907961405680C05680CFF00FFFF00FFFF00FFFF00FFFF00FF
            0A6A150A7F150BB61C09B91A08B41807B21609B31909B41909B81A09B91A0783
            10044D06FF00FFFF00FFFF00FF0B6A150F852216BD3411B7270BB21C07B11608
            B11709B21909B21909B21909B41909BA1A07841006670CFF00FFFF00FF0B6A15
            20BE491BBD4014B7300AB21F28BC36DFF5E1EEFAEF63CE6D09B21909B21909B3
            1909BA1A06670CFF00FF0872101B9A3A2AC65B1DBB450EB4250BB31B11B4219A
            DFA0FFFFFFF7FDF85ACB6509B21909B21909B81A089413045D090872102AB65B
            2CC56522BD4D0FB4220AB21A0CB31C0AB2198DDB95FDFEFDF6FCF758CB6309B2
            1909B51A08AB17045D090F821C37C26C33C76CCDF1DAC9EFD3C7EED0C8EFD2C5
            EED0C7EECFF8FDF9FFFFFFF2FBF36FD27908B41909B31905650B138D2358CC83
            42C977FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDFEFDFFFFFFFFFFFFBCEA
            C10AB41A09B319066D0D0F911D6FD2935FD38D6DD49572D69971D69872D69964
            D28C92DFA8FBFEFBFFFFFFACE5B82EBF4C11B82B08B11905610A0F911D67CC83
            9BE5BA38C67030C36938C56F38C56F70D697E8F8EEFFFFFF9FE2B120BD481AB9
            3E10BA2908A31705610AFF00FF25AE39BCEDD282DBA428C0632FC26753CD82F7
            FDF9FFFFFF9CE2B222BC4B1DBA4118B73614C0300A8517FF00FFFF00FF25AE39
            71D28CD2F4E180DAA336C46D39C56FBCECCEABE6C22DC26324BE5623BC4D1FC1
            4616AE340A8517FF00FFFF00FFFF00FF25AE3984D89FDBF7EAAFE8C66BD49352
            CC8144C97849CA7B48CB7839CB6A21B6490F7C1FFF00FFFF00FFFF00FFFF00FF
            FF00FF25AE3925AE39ADE8C5CCF2DEBAEDD1A6E7C291E2B364D4922FB1572FB1
            57FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF32B74E25AE3925
            AE3925AE3925AE3924A342FF00FFFF00FFFF00FFFF00FFFF00FF}
          ParentFont = False
          TabOrder = 1
          Visible = False
          OnClick = BtnMovementsClick
        end
        object RGKernel: TRadioGroup
          Left = 165
          Top = 24
          Width = 188
          Height = 67
          Caption = 'Dispersal Kernel'
          Ctl3D = True
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ItemIndex = 0
          Items.Strings = (
            'Negative exponential'
            'Double negative exponential')
          ParentCtl3D = False
          ParentFont = False
          TabOrder = 2
          OnClick = RGKernelClick
        end
        object PanelKernels: TPanel
          Left = 9
          Top = 95
          Width = 416
          Height = 302
          Alignment = taLeftJustify
          BevelInner = bvLowered
          Caption = 'Dispersal Kernels'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentFont = False
          TabOrder = 3
          VerticalAlignment = taAlignTop
          DesignSize = (
            416
            302)
          object VarTextPanel: TPanel
            Left = 145
            Top = 141
            Width = 177
            Height = 96
            BevelEdges = [beBottom]
            BevelOuter = bvNone
            TabOrder = 4
            Visible = False
            DesignSize = (
              177
              96)
            object VarLabel1: TLabel
              Left = 33
              Top = 11
              Width = 131
              Height = 16
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              BiDiMode = bdLeftToRight
              Caption = 'Mean distance I (m)'
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentBiDiMode = False
              ParentFont = False
            end
            object VarLabel2: TLabel
              Left = 44
              Top = 41
              Width = 119
              Height = 16
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              BiDiMode = bdLeftToRight
              Caption = 'Mean distance II (m)'
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentBiDiMode = False
              ParentFont = False
            end
            object VarLabel3: TLabel
              Left = 106
              Top = 73
              Width = 54
              Height = 16
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              BiDiMode = bdLeftToRight
              Caption = 'P kernel I'
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentBiDiMode = False
              ParentFont = False
            end
          end
          object edtDist1Scale: TLabeledEdit
            Left = 319
            Top = 151
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.AlignWithMargins = True
            EditLabel.Width = 45
            EditLabel.Height = 32
            EditLabel.BiDiMode = bdLeftToRight
            EditLabel.Caption = 'Scaling factor:'
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentBiDiMode = False
            EditLabel.ParentFont = False
            EditLabel.WordWrap = True
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 5
            Text = '0.0'
            Visible = False
          end
          object edtDist2Scale: TEdit
            Left = 319
            Top = 182
            Width = 50
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 6
            Text = '0.0'
            Visible = False
          end
          object edtPKern1Scale: TEdit
            Left = 319
            Top = 212
            Width = 50
            Height = 24
            Alignment = taRightJustify
            Anchors = []
            AutoSize = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 7
            Text = '0.0'
            Visible = False
          end
          object BtnUpdKernel: TButton
            Left = 70
            Top = 260
            Width = 113
            Height = 25
            Caption = 'Update graph'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ParentFont = False
            TabOrder = 8
            OnClick = BtnUpdKernelClick
          end
          object CBSexKernels: TCheckBox
            Left = 12
            Top = 24
            Width = 133
            Height = 17
            Caption = 'Sex dependent'
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 0
            Visible = False
            OnClick = CBSexKernelsClick
          end
          object CBStageKernels: TCheckBox
            Left = 12
            Top = 48
            Width = 173
            Height = 17
            Caption = 'Stage dependent'
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 1
            Visible = False
            OnClick = CBStageKernelsClick
          end
          object CBIndVarKernel: TCheckBox
            Left = 12
            Top = 72
            Width = 157
            Height = 17
            Caption = 'Individual variability'
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = [fsBold]
            ParentFont = False
            TabOrder = 2
            OnClick = CBIndVarKernelClick
          end
          object PanelDist: TPanel
            Left = -9
            Top = 127
            Width = 322
            Height = 127
            Anchors = []
            BevelEdges = []
            BevelOuter = bvNone
            TabOrder = 3
            DesignSize = (
              322
              127)
            object edtDist1: TLabeledEdit
              Left = 142
              Top = 25
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              EditLabel.Width = 119
              EditLabel.Height = 16
              EditLabel.Caption = 'Mean distance I (m) '
              EditLabel.Font.Charset = DEFAULT_CHARSET
              EditLabel.Font.Color = clWindowText
              EditLabel.Font.Height = -13
              EditLabel.Font.Name = 'Tahoma'
              EditLabel.Font.Style = []
              EditLabel.ParentFont = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              LabelPosition = lpLeft
              ParentFont = False
              TabOrder = 0
              Text = '0.0'
            end
            object edtDist2: TLabeledEdit
              Left = 142
              Top = 56
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              EditLabel.Width = 123
              EditLabel.Height = 16
              EditLabel.Caption = 'Mean distance II (m) '
              EditLabel.Font.Charset = DEFAULT_CHARSET
              EditLabel.Font.Color = clWindowText
              EditLabel.Font.Height = -13
              EditLabel.Font.Name = 'Tahoma'
              EditLabel.Font.Style = []
              EditLabel.ParentFont = False
              Enabled = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              LabelPosition = lpLeft
              ParentFont = False
              TabOrder = 3
              Text = '0.0'
              Visible = False
            end
            object edtPKern1: TLabeledEdit
              Left = 141
              Top = 86
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              EditLabel.Width = 58
              EditLabel.Height = 16
              EditLabel.Caption = 'P kernel I '
              EditLabel.Font.Charset = DEFAULT_CHARSET
              EditLabel.Font.Color = clWindowText
              EditLabel.Font.Height = -13
              EditLabel.Font.Name = 'Tahoma'
              EditLabel.Font.Style = []
              EditLabel.ParentFont = False
              Enabled = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              LabelPosition = lpLeft
              ParentFont = False
              TabOrder = 6
              Text = '0.0'
              Visible = False
            end
            object edtDist2Mean: TEdit
              Left = 198
              Top = 56
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentFont = False
              TabOrder = 4
              Text = '0.0'
              Visible = False
            end
            object edtDist1Mean: TLabeledEdit
              Left = 198
              Top = 25
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              EditLabel.Width = 36
              EditLabel.Height = 16
              EditLabel.Caption = 'Mean:'
              EditLabel.Font.Charset = DEFAULT_CHARSET
              EditLabel.Font.Color = clWindowText
              EditLabel.Font.Height = -13
              EditLabel.Font.Name = 'Tahoma'
              EditLabel.Font.Style = []
              EditLabel.ParentFont = False
              Enabled = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentFont = False
              TabOrder = 1
              Text = '0.0'
              Visible = False
            end
            object edtDist1SD: TLabeledEdit
              Left = 253
              Top = 25
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              EditLabel.Width = 28
              EditLabel.Height = 16
              EditLabel.Caption = 'S.d.:'
              EditLabel.Font.Charset = DEFAULT_CHARSET
              EditLabel.Font.Color = clWindowText
              EditLabel.Font.Height = -13
              EditLabel.Font.Name = 'Tahoma'
              EditLabel.Font.Style = []
              EditLabel.ParentFont = False
              Enabled = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentFont = False
              TabOrder = 2
              Text = '0.0'
              Visible = False
            end
            object edtDist2SD: TEdit
              Left = 253
              Top = 56
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentFont = False
              TabOrder = 5
              Text = '0.0'
              Visible = False
            end
            object edtPKern1SD: TEdit
              Left = 253
              Top = 88
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentFont = False
              TabOrder = 8
              Text = '0.0'
              Visible = False
            end
            object edtPKern1Mean: TEdit
              Left = 197
              Top = 86
              Width = 50
              Height = 24
              Alignment = taRightJustify
              Anchors = []
              AutoSize = False
              Font.Charset = DEFAULT_CHARSET
              Font.Color = clWindowText
              Font.Height = -13
              Font.Name = 'Tahoma'
              Font.Style = []
              ParentFont = False
              TabOrder = 7
              Text = '0.0'
              Visible = False
            end
          end
          object KernelGraph: TChart
            Left = 189
            Top = 48
            Width = 220
            Height = 254
            Legend.Visible = False
            Title.Text.Strings = (
              '')
            Title.Visible = False
            BottomAxis.Automatic = False
            BottomAxis.AutomaticMinimum = False
            BottomAxis.AxisValuesFormat = '###,##0.#'
            BottomAxis.Grid.Visible = False
            BottomAxis.Title.Caption = 'Distance (m)'
            BottomAxis.Title.Font.Style = [fsBold, fsItalic]
            LeftAxis.Automatic = False
            LeftAxis.AutomaticMinimum = False
            LeftAxis.AxisValuesFormat = '0.######'
            LeftAxis.Grid.Visible = False
            LeftAxis.Title.Font.Style = [fsBold, fsItalic]
            RightAxis.Automatic = False
            RightAxis.AutomaticMaximum = False
            RightAxis.AutomaticMinimum = False
            RightAxis.Visible = False
            TopAxis.Automatic = False
            TopAxis.AutomaticMaximum = False
            TopAxis.AutomaticMinimum = False
            TopAxis.Visible = False
            View3D = False
            View3DWalls = False
            TabOrder = 9
            DefaultCanvas = 'TGDIPlusCanvas'
            ColorPaletteIndex = 13
            object Kernel1: TAreaSeries
              SeriesColor = clBlue
              AreaChartBrush.Color = clGray
              AreaChartBrush.BackColor = clDefault
              AreaLinesPen.Visible = False
              DrawArea = True
              Pointer.InflateMargins = True
              Pointer.Style = psRectangle
              Pointer.Visible = False
              XValues.Name = 'X'
              XValues.Order = loAscending
              YValues.Name = 'Y'
              YValues.Order = loNone
            end
            object Kernel2: TAreaSeries
              SeriesColor = clRed
              AreaChartBrush.Color = clGray
              AreaChartBrush.BackColor = clDefault
              AreaLinesPen.Visible = False
              DrawArea = True
              Pointer.InflateMargins = True
              Pointer.Style = psRectangle
              Pointer.Visible = False
              Transparency = 50
              XValues.Name = 'X'
              XValues.Order = loAscending
              YValues.Name = 'Y'
              YValues.Order = loNone
            end
          end
        end
        object PanelMortality: TPanel
          Left = 9
          Top = 403
          Width = 416
          Height = 118
          Alignment = taLeftJustify
          BevelInner = bvLowered
          Caption = 'Dispersal Mortality'
          TabOrder = 4
          VerticalAlignment = taAlignTop
          object edtMortProb: TLabeledEdit
            Left = 319
            Top = 21
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 116
            EditLabel.Height = 16
            EditLabel.Caption = 'Mortality probability '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 1
            Text = '0.0'
          end
          object edtMortSlope: TLabeledEdit
            Left = 320
            Top = 51
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 38
            EditLabel.Height = 16
            EditLabel.Caption = 'slope  '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 2
            Text = '0.0'
            Visible = False
          end
          object edtMortInfl: TLabeledEdit
            Left = 320
            Top = 81
            Width = 50
            Height = 24
            Alignment = taRightJustify
            EditLabel.Width = 87
            EditLabel.Height = 16
            EditLabel.Caption = 'inflection point '
            EditLabel.Font.Charset = DEFAULT_CHARSET
            EditLabel.Font.Color = clWindowText
            EditLabel.Font.Height = -13
            EditLabel.Font.Name = 'Tahoma'
            EditLabel.Font.Style = []
            EditLabel.ParentFont = False
            Enabled = False
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            LabelPosition = lpLeft
            ParentFont = False
            TabOrder = 3
            Text = '0.0'
            Visible = False
          end
          object RGMortality: TRadioGroup
            Left = 12
            Top = 24
            Width = 150
            Height = 73
            Ctl3D = True
            Font.Charset = DEFAULT_CHARSET
            Font.Color = clWindowText
            Font.Height = -13
            Font.Name = 'Tahoma'
            Font.Style = []
            ItemIndex = 0
            Items.Strings = (
              'Constant'
              'Distance Dependent')
            ParentCtl3D = False
            ParentFont = False
            TabOrder = 0
            OnClick = RGMortalityClick
          end
        end
      end
    end
    object TSDispersal2: TTabSheet
      Caption = 'Sex-/stage-dependent dispersal'
      ImageIndex = 2
      OnEnter = TSDispersal2Enter
      ExplicitLeft = 0
      ExplicitTop = 0
      ExplicitWidth = 0
      ExplicitHeight = 0
      object PanelSexEmig: TPanel
        Left = 0
        Top = 5
        Width = 897
        Height = 174
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 0
        VerticalAlignment = taAlignTop
        Visible = False
        object LabelSexEmigPar: TLabel
          Left = 9
          Top = 18
          Width = 137
          Height = 16
          AutoSize = False
          Caption = 'Density independent:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          WordWrap = True
        end
        object LabelSexEmigDensPar: TLabel
          Left = 263
          Top = 18
          Width = 137
          Height = 16
          AutoSize = False
          Caption = 'Density dependent:'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          Visible = False
          WordWrap = True
        end
        object SexEmigLabel: TLabel
          Left = 8
          Top = 2
          Width = 87
          Height = 16
          Caption = 'EMIGRATION  '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object SexEmigPar: TStringGrid
          Left = 9
          Top = 40
          Width = 239
          Height = 129
          BiDiMode = bdLeftToRight
          ColCount = 3
          DefaultColWidth = 70
          DoubleBuffered = False
          RowCount = 3
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 0
          ColWidths = (
            70
            70
            70)
          RowHeights = (
            24
            24
            24)
        end
        object SexEmigDensPar: TStringGrid
          Left = 263
          Top = 40
          Width = 625
          Height = 129
          BiDiMode = bdLeftToRight
          ColCount = 4
          DefaultColWidth = 70
          DoubleBuffered = False
          Enabled = False
          RowCount = 3
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 1
          Visible = False
          ColWidths = (
            70
            70
            70
            70)
          RowHeights = (
            24
            24
            24)
        end
      end
      object PanelSexTransfer: TPanel
        Left = 0
        Top = 183
        Width = 897
        Height = 174
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 1
        VerticalAlignment = taAlignTop
        Visible = False
        object SexTransferLabel: TLabel
          Left = 8
          Top = 2
          Width = 205
          Height = 16
          Caption = 'TRANSFER - Dispersal Kernels    '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object SexKernPar: TStringGrid
          Left = 9
          Top = 21
          Width = 880
          Height = 145
          BiDiMode = bdLeftToRight
          DefaultColWidth = 100
          DoubleBuffered = False
          RowCount = 3
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 0
          ColWidths = (
            100
            100
            100
            100
            100)
          RowHeights = (
            24
            24
            24)
        end
      end
      object PanelSexSettle: TPanel
        Left = 0
        Top = 363
        Width = 897
        Height = 174
        Alignment = taLeftJustify
        BevelInner = bvLowered
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = [fsBold, fsItalic]
        ParentFont = False
        TabOrder = 2
        VerticalAlignment = taAlignTop
        Visible = False
        object SexSettleParLabel: TLabel
          Left = 8
          Top = 2
          Width = 240
          Height = 16
          Caption = 'SETTLEMENT - Movement Processes   '
          Color = clSkyBlue
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clBlue
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = [fsBold, fsItalic]
          ParentColor = False
          ParentFont = False
        end
        object SexSettlePar: TStringGrid
          Left = 8
          Top = 24
          Width = 880
          Height = 145
          BiDiMode = bdLeftToRight
          ColCount = 9
          DefaultColWidth = 95
          DoubleBuffered = False
          RowCount = 3
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -13
          Font.Name = 'Tahoma'
          Font.Style = []
          Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
          ParentBiDiMode = False
          ParentDoubleBuffered = False
          ParentFont = False
          TabOrder = 0
          ColWidths = (
            95
            95
            95
            95
            95
            95
            95
            95
            95)
          RowHeights = (
            24
            24
            24)
        end
      end
    end
  end
  object BtnOK: TButton
    Left = 812
    Top = 575
    Width = 89
    Height = 25
    Caption = 'OK'
    Default = True
    TabOrder = 1
    OnClick = BtnOKClick
  end
  object BtnCancel: TButton
    Left = 907
    Top = 575
    Width = 89
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    ModalResult = 2
    TabOrder = 2
    OnClick = BtnCancelClick
  end
end
