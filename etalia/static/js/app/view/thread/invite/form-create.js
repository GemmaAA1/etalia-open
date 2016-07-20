define(['app', 'app/model/thread/invite', 'select2'], function (App) {

    App.View.Thread.InviteCreateForm = App.Backbone.Form.extend({

        schema: {
            to_user: {
                type: 'Select',
                title: 'Who',
                options: [],
                validators: ['required', App.Model.Invite.validators.to_user]
            }
        },

        template: App._.template('\
            <form class="thread-create-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <div class="form-group form-footer buttons">\
                    <div class="col-sm-offset-2 col-sm-10">\
                        <button type="submit" class="btn btn-primary form-submit">\
                            Invite\
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
        },

        postRender: function() {

            function formatResult(state) {
                //console.log(state);
                if (!state.id) {
                    return $('<span>' + state.text + '</span>');
                }

                return $('<span>' + state.first_name + ' ' + state.last_name + '</span>');
            }

            App.$('select[name="to_user"]').select2({
                ajax: {
                    url: App.config.api_root + "/user/users/",
                    dataType: 'json',
                    delay: 250,
                    data: function (params) {
                        return {
                            search: params.term, // search term
                            page: params.page
                        };
                    },
                    processResults: function (data, params) {
                        // parse the results into the format expected by Select2
                        // since we are using custom formatting functions we do not need to
                        // alter the remote JSON data, except to indicate that infinite
                        // scrolling can be used
                        params.page = params.page || 1;

                        return {
                            results: data.results,
                            pagination: {
                                more: (params.page * 30) < data.total_count
                            }
                        };
                    },
                    cache: true
                },
                escapeMarkup: function (markup) { return markup; }, // let our custom formatter work
                minimumInputLength: 3,
                templateResult: formatResult, // omitted for brevity, see the source of this page
                templateSelection: formatResult // omitted for brevity, see the source of this page
            });

            return this;
        },

        _onSubmit: function(e) {
            e.preventDefault();

            var errors = this.commit();
            if (errors) {
                App.log('Errors', errors);

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


    App.View.Thread.InviteCreateForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = App.Model.Invite.createNew();
        }

        var form = new App.View.Thread.InviteCreateForm(options);
        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.InviteCreateForm;
});
