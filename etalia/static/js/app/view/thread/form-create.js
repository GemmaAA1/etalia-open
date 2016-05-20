define(['app', 'app/model/thread/thread'], function (App) {

    App.View.Thread = App.View.Thread || {};

    App.View.Thread.CreateForm = App.Backbone.Form.extend({

        schema: {
            privacy: {
                type: "Radio",
                options: [
                    {label: "Public", val: App.Model.Thread.PRIVACY_PUBLIC},
                    {label: "Private", val: App.Model.Thread.PRIVACY_PRIVATE}
                ],
                help: 'Interdum et malesuada fames ac ante ipsum primis in faucibus. ' +
                    'Sed volutpat ante ut sodales pellentesque. Sed at est sed diam tempus molestie. ' +
                    'Integer sit amet egestas tortor.',
                validators: ['required', App.Model.Thread.validators.privacy]
            },

            type: {
                type: "Radio",
                options: [
                    { label: "Question", val: App.Model.Thread.TYPE_QUESTION},
                    { label: "Paper", val: App.Model.Thread.TYPE_PAPER}
                ],
                validators: ['required', App.Model.Thread.validators.type]
            },

            paper: {
                type: 'Select',
                options: function(callback) {
                    var user = App.getCurrentUser(),
                        papers = new App.Model.Papers();

                    papers.url = App.config.api_root + '/user/user-libs/' + user.get('id') + '/papers';
                    papers
                        .fetch()
                        .then(function() {
                            var options = [{val: null, label: 'None'}];
                            papers.each(function (paper) {
                                options.push({val: paper.get('id'), label: paper.get('title')});
                            });
                            callback(options);
                        }, function(jqXHR, textStatus, errorThrown) {
                            throw textStatus + ' ' + errorThrown;
                        });

                },
                validators: [App.Model.Thread.validators.paper]
            },

            title: {
                type: 'Text',
                validators: ['required', App.Model.Thread.validators.title]
            }
        },

        template: App._.template('\
            <form class="thread-create-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <div class="form-group form-footer buttons">\
                    <div class="col-sm-offset-2 col-sm-10">\
                        <button type="submit" class="btn btn-primary form-submit">\
                            Create\
                        </button>\
                        <button type="button" class="btn btn-default form-cancel">\
                            Cancel\
                        </button>\
                    </div>\
                </div>\
            </form>\
        '),

        events: {
            "click .form-submit": "_onSubmit",
            "click .form-cancel": "_onCancel"
        },

        initialize: function () {
            App.Backbone.Form.prototype.initialize.apply(this, arguments);

            this.listenTo(this, "submit", this._onSubmit);

            this.listenTo(this, "type:change", this._onTypeChange);
        },

        render: function () {
            App.Backbone.Form.prototype.render.apply(this, arguments);

            if (this.getEditor('type').getValue() == App.Model.Thread.TYPE_PAPER) {
                this.$('.field-paper').show();
            } else {
                this.$('.field-paper').hide();
            }

            return this;
        },

        _onTypeChange: function(form, typeEditor) {
            if (typeEditor.getValue() == App.Model.Thread.TYPE_PAPER) {
                this.$('.field-paper').show();
            } else {
                this.$('.field-paper').hide();
            }
        },

        _onSubmit: function(e) {
            e.preventDefault();

            var errors = this.commit();
            if (errors) {
                App.log('validation errors', errors);

                this.trigger('validation_failure', this);
                return;
            }

            this.trigger('validation_success', this);
        },

        _onCancel: function(e) {
            e.preventDefault();

            this.trigger('cancel', this);
        }
    });


    App.View.Thread.CreateForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = App.Model.Thread.createNew();
        }

        var form = new App.View.Thread.CreateForm(options);
        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.CreateForm;
});
