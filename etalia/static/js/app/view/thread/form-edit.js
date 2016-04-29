define([
    'app',
    'app/model/thread/thread'
], function (App) {

    App.View.Thread = App.View.Thread || {};

    App.View.Thread.EditForm = App.Backbone.Form.extend({

        schema: function() {
            var result = {
                privacy: {
                    type: "Radio",
                    options: [
                        { label: "Public", val: App.Model.Thread.PRIVACY_PUBLIC},
                        { label: "Private", val: App.Model.Thread.PRIVACY_PRIVATE}
                    ],
                    help: 'Interdum et malesuada fames ac ante ipsum primis in faucibus. ' +
                    'Sed volutpat ante ut sodales pellentesque. Sed at est sed diam tempus molestie. ' +
                    'Integer sit amet egestas tortor.'
                }
            };

            if (this.model.get('type') === App.Model.Thread.TYPE_PAPER) {
                result.paper = {
                    type: 'Select',
                    options: function (callback) {
                        var user = App.getCurrentUser(),
                            papers = new App.Collection.Papers();

                        papers.url = App.config.api_root + '/user/user-libs/' + user.get('id') + '/papers';
                        papers
                            .fetch()
                            .then(function () {
                                var options = [{val: null, label: 'None'}];
                                papers.each(function (paper) {
                                    options.push({val: paper.get('id'), label: paper.get('title')});
                                });
                                callback(options);
                            }, function (jqXHR, textStatus, errorThrown) {
                                throw textStatus + ' ' + errorThrown;
                            });
                    }
                    // TODO validation
                }
            }

            result.title = {type: 'Text'};

            return result;
        },

        template: App._.template('\
            <form class="thread-create-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <div class="form-group form-footer">\
                    <div class="col-sm-offset-2 col-sm-10">\
                        <button type="submit" class="btn btn-primary form-submit">\
                            Update\
                        </button>\
                        <button type="button" class="btn btn-default form-cancel">\
                            Cancel\
                        </button>\
                    </div>\
                </div>\
            </form>\
        '),

        events: {
            "click .form-submit": "onSubmit",
            "click .form-cancel": "onCancel"
        },

        initialize: function () {
            App.Backbone.Form.prototype.initialize.apply(this, arguments);

            this.listenTo(this, "submit", this.onSubmit);
        },

        onSubmit: function(e) {
            e.preventDefault();

            var errors = this.commit({ validate: true });
            if (errors) {
                App.log('Errors', errors);

                this.trigger('validation_failure', this);
                return;
            }

            this.trigger('validation_success', this);
        },

        onCancel: function(e) {
            e.preventDefault();

            this.trigger('cancel', this);
        }
    });


    App.View.Thread.EditForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'Expected options.model instance of App.Model.Thread';
        }

        var form = new App.View.Thread.EditForm(options);
        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.EditForm;
});